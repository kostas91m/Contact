using System;
using System.Collections.Generic;
using MGroup.Constitutive.Structural;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.DataStructures;

namespace MGroup.FEM.Structural.Line
{
	public class ContactNodeToSurface3DFriction : IStructuralElementType
	{
		private static readonly IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ };
		private static readonly IDofType[][] dofs = new IDofType[][] { nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes };
		private readonly double PenaltyFactorNormal;
		private readonly double PenaltyFactorTangential;
		private readonly double StickingCoefficient;
		private readonly double SlidingCoefficient;
		private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
		private double[] DisplacementVector { get; set; }
		private double ContactArea { get; }
		private Dictionary<int, double[]> StickingPoint { get; set; }
		private Dictionary<int, double[]> TangentialTraction { get; set; }
		private Dictionary<int, double[]> CPPSurfaceBaseVector1 { get; set; }
		private Dictionary<int, double[]> CPPSurfaceBaseVector2 { get; set; }
		private double[] PreviousConvergedSolutionNodalCoordinates { get; set; }
		public ContactNodeToSurface3DFriction(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplierNormal,
			double penaltyFactorMultiplierTangential, double stickingCoefficient, double slidingCoefficient,
			double contactArea)
		{
			if (nodes.Count != 5)
			{
				throw new ArgumentException("This Constructor can only be used for linear Surfaces");
			}
			this.PenaltyFactorNormal = penaltyFactorMultiplierNormal * youngModulus;
			this.PenaltyFactorTangential = penaltyFactorMultiplierTangential * youngModulus;
			this.StickingCoefficient = stickingCoefficient;
			this.SlidingCoefficient = slidingCoefficient;
			this.DisplacementVector = new double[15];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
			this.StickingPoint = new Dictionary<int, double[]>()
			{
				{ 1, new double[2] },
			};
			this.TangentialTraction = new Dictionary<int, double[]>()
			{
				{ 1, new double[2] },
			};
			this.CPPSurfaceBaseVector1 = new Dictionary<int, double[]>()
			{
				{ 1, new double[3] },
			};
			this.CPPSurfaceBaseVector2 = new Dictionary<int, double[]>()
			{
				{ 1, new double[3] },
			};
			this.PreviousConvergedSolutionNodalCoordinates = new double[3 * nodes.Count];
			for (var i = 0; i < nodes.Count; i++)
			{
				this.PreviousConvergedSolutionNodalCoordinates[3 * i] = nodes[i].X;
				this.PreviousConvergedSolutionNodalCoordinates[3 * i + 1] = nodes[i].Y;
				this.PreviousConvergedSolutionNodalCoordinates[3 * i + 2] = nodes[i].Z;
			}
			InitializeTangentialProperties();
		}
		public ContactNodeToSurface3DFriction(IReadOnlyList<INode> nodes, double penaltyFactorNormal,
			double penaltyFactorTangential, double stickingCoefficient, double slidingCoefficient,
			double contactArea)
		{
			if (nodes.Count != 5)
			{
				throw new ArgumentException("This Constructor can only be used for linear Surfaces");
			}
			this.PenaltyFactorNormal = penaltyFactorNormal;
			this.PenaltyFactorTangential = penaltyFactorTangential;
			this.StickingCoefficient = stickingCoefficient;
			this.SlidingCoefficient = slidingCoefficient;
			this.DisplacementVector = new double[15];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
			this.StickingPoint = new Dictionary<int, double[]>()
			{
				{ 1, new double[2] },
			};
			this.TangentialTraction = new Dictionary<int, double[]>()
			{
				{ 1, new double[2] },
			};
			this.CPPSurfaceBaseVector1 = new Dictionary<int, double[]>()
			{
				{ 1, new double[3] },
			};
			this.CPPSurfaceBaseVector2 = new Dictionary<int, double[]>()
			{
				{ 1, new double[3] },
			};
			this.PreviousConvergedSolutionNodalCoordinates = new double[3 * nodes.Count];
			for (var i = 0; i < nodes.Count; i++)
			{
				this.PreviousConvergedSolutionNodalCoordinates[3 * i] = nodes[i].X;
				this.PreviousConvergedSolutionNodalCoordinates[3 * i + 1] = nodes[i].Y;
				this.PreviousConvergedSolutionNodalCoordinates[3 * i + 2] = nodes[i].Z;
			}
			InitializeTangentialProperties();
		}
		public ContactNodeToSurface3DFriction(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplierNormal,
			double penaltyFactorMultiplierTangential, double stickingCoefficient, double slidingCoefficient,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, youngModulus, penaltyFactorMultiplierNormal, penaltyFactorMultiplierTangential,
				   stickingCoefficient, slidingCoefficient, 1d)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public ContactNodeToSurface3DFriction(IReadOnlyList<INode> nodes, double penaltyFactorNormal,
			double penaltyFactorTangential, double stickingCoefficient, double slidingCoefficient,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, penaltyFactorNormal, penaltyFactorTangential,
				   stickingCoefficient, slidingCoefficient, 1d)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public CellType CellType { get; } = CellType.Unknown;

		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		private double[] XUpdatedVector()
		{
			var xVectorUpdated = new double[15];
			for (var i = 0; i < 5; i++)
			{
				xVectorUpdated[3 * i] = Nodes[i].X + DisplacementVector[3 * i];
				xVectorUpdated[3 * i + 1] = Nodes[i].Y + DisplacementVector[3 * i + 1];
				xVectorUpdated[3 * i + 2] = Nodes[i].Z + DisplacementVector[3 * i + 2];
			}
			return xVectorUpdated;
		}

		private Dictionary<int, double[]> SurfaceVectors(double[,] da1, double[,] da2)
		{
			var xUpdated = XUpdatedVector();
			return new Dictionary<int, double[]>()
			{
				{ 1, Matrix.CreateFromArray(da1).Multiply(xUpdated).Scale(-1.0) },
				{ 2, Matrix.CreateFromArray(da2).Multiply(xUpdated).Scale(-1.0) }
			};
		}

		private double[,] MetricTensor(Dictionary<int, double[]> masterSurfaceVectors) => new double[,]
		{
			{ masterSurfaceVectors[1].DotProduct(masterSurfaceVectors[1]), masterSurfaceVectors[1].DotProduct(masterSurfaceVectors[2]) },
			{ masterSurfaceVectors[2].DotProduct(masterSurfaceVectors[1]), masterSurfaceVectors[2].DotProduct(masterSurfaceVectors[2]) }
		};

		private double MetricTensorDet(double[,] m)
		{
			var detm = m[0, 0] * m[1, 1] - m[1, 0] * m[0, 1];
			return detm;
		}

		private Matrix InverseMetricTensor(double[,] m)
		{
			var detm = MetricTensorDet(m);
			var mInv = Matrix.CreateFromArray(new double[,] { { m[1, 1], -m[0, 1] }, { -m[1, 0], m[0, 0] } }).Scale(1.0 / detm);
			return mInv;
		}

		private double[] NormalVector(double[,] metricTensor, Dictionary<int, double[]> masterSurfaceVectors)
		{
			var n = (masterSurfaceVectors[1].CrossProduct(masterSurfaceVectors[2])).Scale(1.0 / (Math.Sqrt(MetricTensorDet(metricTensor))));
			return n;
		}

		private double[] Calculate_f(Dictionary<int, double[]> masterSurfaceVectors, double[,] aMatrix, double[] xUpdated) => new double[]
		{
			Vector.CreateFromArray(masterSurfaceVectors[1]).DotProduct(Matrix.CreateFromArray(aMatrix) * Vector.CreateFromArray(xUpdated)),
			Vector.CreateFromArray(masterSurfaceVectors[2]).DotProduct(Matrix.CreateFromArray(aMatrix) * Vector.CreateFromArray(xUpdated))
		};

		private double Calculate_e(double[,] aMatrix, double[] xUpdated)
		{
			var e = Vector.CreateFromArray(new double[]
			{
				0.25*(xUpdated[0] - xUpdated[3] + xUpdated[6] - xUpdated[9]),
				0.25*(xUpdated[1] - xUpdated[4] + xUpdated[7] - xUpdated[10]),
				0.25*(xUpdated[2] - xUpdated[5] + xUpdated[8] - xUpdated[11])
			}).DotProduct(Matrix.CreateFromArray(aMatrix) * Vector.CreateFromArray(xUpdated));
			return e;
		}

		private Vector CalculateDeltaKsi(double detm, double[,] mTensor, double[] fVector, double e)
		{
			var scalar = 1.0 / (detm - Math.Pow(e, 2) + 2.0 * e * mTensor[0, 1]);
			var matrix = new double[,]
			{
				{mTensor[1,1], e-mTensor[0,1] },
				{e-mTensor[1,0], mTensor[0,0] }
			};
			var deltaKsi = (Matrix.CreateFromArray(matrix) * Vector.CreateFromArray(fVector)).Scale(scalar);
			return deltaKsi;
		}

		private Vector Project(double[] ksiVectorInitial)
		{
			var maxIterations = 1000;
			var tol = Math.Pow(10.0, -4.0);
			var norm = new double();
			var ksiVector = Vector.CreateFromArray(ksiVectorInitial);
			var xUpdated = XUpdatedVector();
			for (var i = 1; i <= maxIterations; i++)
			{
				var aMatrices = CalculatePositionMatrix(ksiVector[0], ksiVector[1]);
				var masterSurfaceVectors = SurfaceVectors(aMatrices.Item2, aMatrices.Item3);
				var f = Calculate_f(masterSurfaceVectors, aMatrices.Item1, xUpdated);
				var e = Calculate_e(aMatrices.Item1, xUpdated);
				var m = MetricTensor(masterSurfaceVectors);
				var detm = MetricTensorDet(m);
				var deltaKsi = CalculateDeltaKsi(detm, m, f, e);
				ksiVector += deltaKsi;
				norm = (deltaKsi).Norm2();
				if (norm <= tol)
				{
					break;
				}
			}
			if (norm > tol)
			{
				throw new Exception("CPP not found in current iterations");
			}
			else
			{
				return ksiVector;
			}

		}

		private Tuple<double[,], double[,], double[,]> CalculatePositionMatrix(double ksi1, double ksi2)
		{
			var N1 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 - ksi2);
			var N2 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 - ksi2);
			var N3 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 + ksi2);
			var N4 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 + ksi2);

			var dN11 = -1.0 / 4.0 * (1.0 - ksi2);
			var dN21 = 1.0 / 4.0 * (1.0 - ksi2);
			var dN31 = 1.0 / 4.0 * (1.0 + ksi2);
			var dN41 = -1.0 / 4.0 * (1.0 + ksi2);

			var dN12 = -1.0 / 4.0 * (1.0 - ksi1);
			var dN22 = -1.0 / 4.0 * (1.0 + ksi1);
			var dN32 = 1.0 / 4.0 * (1.0 + ksi1);
			var dN42 = 1.0 / 4.0 * (1.0 - ksi1);

			var aMatrix = new double[,]
				{
					{ -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, 1.0, 0.0, 0.0 },
					{ 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, 1.0, 0.0 },
					{ 0.0, 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, 1.0 }
				};

			var da1Matrix = new double[,]
				{
					{ -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0 }
				};

			var da2Matrix = new double[,]
				{
					{ -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0 }
				};
			return new Tuple<double[,], double[,], double[,]>(aMatrix, da1Matrix, da2Matrix);
		}

		private double CalculatePenetration(double[] normalVector, double[,] aMatrix, double[] xUpdated)
		{
			var ksi3 = Vector.CreateFromArray(xUpdated).DotProduct(Matrix.CreateFromArray(aMatrix).Transpose() * Vector.CreateFromArray(normalVector));
			return ksi3;
		}

		private Tuple<double[], double[]> MasterPrevSurfaceBaseVectors(double[,] da1Matrix, double[,] da2Matrix, double[] xupd)
		{
			var surfaceVector1 = Matrix.CreateFromArray(da1Matrix).Multiply(xupd);
			var surfaceVector2 = Matrix.CreateFromArray(da2Matrix).Multiply(xupd);
			return new Tuple<double[], double[]>(surfaceVector1, surfaceVector2);
		}

		public void InitializeTangentialProperties()
		{
			var x0 = XUpdatedVector().Scale(-1d);
			var cppInitial = Project(new double[2] { 0.0, 0.0 });
			StickingPoint[1][0] = cppInitial[0];
			StickingPoint[1][1] = cppInitial[1];
			var positionMatrices = CalculatePositionMatrix(cppInitial[0], cppInitial[1]);
			var da1Matrix = positionMatrices.Item2;
			var da2Matrix = positionMatrices.Item3;
			var rM = MasterPrevSurfaceBaseVectors(da1Matrix, da2Matrix, x0);
			CPPSurfaceBaseVector1[1][0] = rM.Item1[0];
			CPPSurfaceBaseVector1[1][1] = rM.Item1[1];
			CPPSurfaceBaseVector1[1][2] = rM.Item1[2];
			CPPSurfaceBaseVector2[1][0] = rM.Item2[0];
			CPPSurfaceBaseVector2[1][1] = rM.Item2[1];
			CPPSurfaceBaseVector2[1][2] = rM.Item2[2];
		}

		public void UpdateTangentialProperties()
		{
			var xUpd = XUpdatedVector();
			var cPP = Project(new double[2] { 0.0, 0.0 });
			if (Math.Abs(cPP[0]) <= 1.05 && Math.Abs(cPP[1]) <= 1.05)
			{
				var positionMatrices = CalculatePositionMatrix(cPP[0], cPP[1]);
				var aMatrix = positionMatrices.Item1;
				var masterSurfaceVectors = SurfaceVectors(positionMatrices.Item2, positionMatrices.Item3);
				var m = MetricTensor(masterSurfaceVectors);
				var n = NormalVector(m, masterSurfaceVectors); 
				var mInv = InverseMetricTensor(m);
				var dRho1 = masterSurfaceVectors[1];
				var dRho2 = masterSurfaceVectors[2];
				var ksi3 = CalculatePenetration(n, positionMatrices.Item1, xUpd);
				if (ksi3 <= 0)
				{
					var xM = XUpdatedVector().Scale(-1d);
					var aMatrixOld = CalculatePositionMatrix(StickingPoint[1][0], StickingPoint[1][1]).Item1;
					var positionVector = new double[3];
					for (var k = 0; k < 4; k++)
					{
						positionVector[0] += xM[3 * k] * aMatrix[0, 3 * k];
						positionVector[1] += xM[3 * k + 1] * aMatrix[0, 3 * k];
						positionVector[2] += xM[3 * k + 2] * aMatrix[0, 3 * k];
					}
					var oldPositionVector = new double[3];
					var nodalDIsplacementsIncrements = CalculateNodalDisplacementsIncrements(xUpd);
					for (var k = 0; k < 4; k++)
					{
						oldPositionVector[0] += -PreviousConvergedSolutionNodalCoordinates[3 * k] * aMatrixOld[0, 3 * k] - nodalDIsplacementsIncrements[3 * k] * aMatrixOld[0, 3 * k];
						oldPositionVector[1] += -PreviousConvergedSolutionNodalCoordinates[3 * k + 1] * aMatrixOld[0, 3 * k] - nodalDIsplacementsIncrements[3 * k + 1] * aMatrixOld[0, 3 * k];
						oldPositionVector[2] += -PreviousConvergedSolutionNodalCoordinates[3 * k + 2] * aMatrixOld[0, 3 * k] - nodalDIsplacementsIncrements[3 * k + 2] * aMatrixOld[0, 3 * k];
					}
					var deltaKsi = new double[]
					{
								(positionVector.Subtract(oldPositionVector)).DotProduct(dRho1) * mInv[0,0] +
								(positionVector.Subtract(oldPositionVector)).DotProduct(dRho2) * mInv[0,1],
								(positionVector.Subtract(oldPositionVector)).DotProduct(dRho1) * mInv[1,0] +
								(positionVector.Subtract(oldPositionVector)).DotProduct(dRho2) * mInv[1,1]
					};
					var rM1 = new double[]
					{
								CPPSurfaceBaseVector1[1][0],
								CPPSurfaceBaseVector1[1][1],
								CPPSurfaceBaseVector1[1][2]
					};
					var rM2 = new double[]
					{
								CPPSurfaceBaseVector2[1][0],
								CPPSurfaceBaseVector2[1][1],
								CPPSurfaceBaseVector2[1][2]
					};
					var mPrev = new double[,]
					{
								{ rM1.DotProduct(rM1), rM1.DotProduct(rM2) },
								{ rM2.DotProduct(rM1), rM2.DotProduct(rM2) }
					};
					var detmPrev = mPrev[0, 0] * mPrev[1, 1] - mPrev[0, 1] * mPrev[1, 0];
					var mPrevInv = new double[,]
					{
								{ mPrev[1,1]/detmPrev, - mPrev[0,1]/detmPrev },
								{ - mPrev[1,0]/detmPrev, mPrev[0,0]/detmPrev }
					};
					var trialTangentialTraction = new double[]
					{
								TangentialTraction[1][0] * mPrevInv[0, 0] * rM1.DotProduct(dRho1) +
								TangentialTraction[1][0] * mPrevInv[0, 1] * rM2.DotProduct(dRho1) +
								TangentialTraction[1][1] * mPrevInv[1, 0] * rM1.DotProduct(dRho1) +
								TangentialTraction[1][1] * mPrevInv[1, 1] * rM2.DotProduct(dRho1) -
								mInv[0, 0] * PenaltyFactorTangential * deltaKsi[0] - mInv[0, 1] * PenaltyFactorTangential * deltaKsi[1],
								TangentialTraction[1][0] * mPrevInv[0, 0] * rM1.DotProduct(dRho2) +
								TangentialTraction[1][0] * mPrevInv[0, 1] * rM2.DotProduct(dRho2) +
								TangentialTraction[1][1] * mPrevInv[1, 0] * rM1.DotProduct(dRho2) +
								TangentialTraction[1][1] * mPrevInv[1, 1] * rM2.DotProduct(dRho2) +
								-mInv[1, 0] * PenaltyFactorTangential * deltaKsi[0] - mInv[1, 1] * PenaltyFactorTangential * deltaKsi[1]
					};
					var tangentialTractionNorm = Math.Pow(mInv[0, 0] * trialTangentialTraction[0] * trialTangentialTraction[0] +
						mInv[0, 1] * trialTangentialTraction[0] * trialTangentialTraction[1] +
						mInv[1, 0] * trialTangentialTraction[1] * trialTangentialTraction[0] +
						mInv[1, 1] * trialTangentialTraction[1] * trialTangentialTraction[1],
						0.5);
					var phiTr = tangentialTractionNorm - StickingCoefficient * PenaltyFactorNormal * Math.Abs(ksi3);
					if (phiTr <= 0)
					{
						TangentialTraction[1][0] = trialTangentialTraction[0];
						TangentialTraction[1][1] = trialTangentialTraction[1];
						StickingPoint[1][0] = cPP[0];
						StickingPoint[1][1] = cPP[1];
					}
					else
					{
						TangentialTraction[1][0] = (trialTangentialTraction[0] / tangentialTractionNorm) * SlidingCoefficient * PenaltyFactorNormal * Math.Abs(ksi3);
						TangentialTraction[1][1] = (trialTangentialTraction[1] / tangentialTractionNorm) * SlidingCoefficient * PenaltyFactorNormal * Math.Abs(ksi3);
						StickingPoint[1][0] = cPP[0];
						StickingPoint[1][1] = cPP[1];
					};
				}
				else
				{
					TangentialTraction[1] = new double[2];
					StickingPoint[1][0] = cPP[0];
					StickingPoint[1][1] = cPP[1];
				}
			}
			else
			{
				TangentialTraction[1] = new double[2];
				StickingPoint[1][0] = cPP[0];
				StickingPoint[1][1] = cPP[1];
			}
			for (var i = 0; i < Nodes.Count; i++)
			{
				PreviousConvergedSolutionNodalCoordinates[i * 3] = xUpd[i * 3];
				PreviousConvergedSolutionNodalCoordinates[i * 3 + 1] = xUpd[i * 3 + 1];
				PreviousConvergedSolutionNodalCoordinates[i * 3 + 2] = xUpd[i * 3 + 2];
			}
			UpdateMasterSurfaceBaseVectors();
		}
		private double[] CalculateNodalDisplacementsIncrements(double[] xUpdated)
		{
			var uNodalIncrements = xUpdated.Subtract(PreviousConvergedSolutionNodalCoordinates);
			return uNodalIncrements;
		}
		private void UpdateMasterSurfaceBaseVectors()
		{
			var ksi = Project(new double[2] { 0.0, 0.0 });
			var positionMatrices = CalculatePositionMatrix(ksi[0], ksi[1]);
			var masterSurfaceVectors = SurfaceVectors(positionMatrices.Item2, positionMatrices.Item3);
			var dRho1 = masterSurfaceVectors[1];
			var dRho2 = masterSurfaceVectors[2];
			CPPSurfaceBaseVector1[1][0] = dRho1[0];
			CPPSurfaceBaseVector1[1][1] = dRho1[1];
			CPPSurfaceBaseVector1[1][2] = dRho1[2];
			CPPSurfaceBaseVector2[1][0] = dRho2[0];
			CPPSurfaceBaseVector2[1][1] = dRho2[1];
			CPPSurfaceBaseVector2[1][2] = dRho2[2];
		}
		#region IElementType Members

		public int ID { get; set; }
		public int SubdomainID { get; set; }
		public IReadOnlyList<INode> Nodes { get; }

		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofs;

		public IReadOnlyList<INode> GetNodesForMatrixAssembly() => Nodes;

		private Matrix CalculateMainStiffnessPart(double[] n, Matrix aMatrix)
		{
			var nxn = n.TensorProduct(n);
			var nxn_A = nxn.MultiplyRight(aMatrix);
			var AT_nxn_A = aMatrix.Transpose().MultiplyRight(nxn_A);
			var mainStiffnessMatrix = AT_nxn_A.Scale(PenaltyFactorNormal);
			return mainStiffnessMatrix;
		}

		private Matrix CalculateRotationalStiffnessPart(double[] normalVector, Matrix aMatrix, Matrix dAd1Matrix, Matrix dAd2Matrix, Dictionary<int, double[]> dRho, double ksi3, double[,] metricTensor)
		{
			var mInv = InverseMetricTensor(metricTensor);
			var scalar1 = PenaltyFactorNormal * ksi3 * mInv[0, 0];
			var scalar2 = PenaltyFactorNormal * ksi3 * mInv[1, 0];
			var scalar3 = PenaltyFactorNormal * ksi3 * mInv[0, 1];
			var scalar4 = PenaltyFactorNormal * ksi3 * mInv[1, 1];
			var mat1 = dAd1Matrix.Transpose() * (normalVector.TensorProduct(dRho[1])) * aMatrix + aMatrix.Transpose() * (dRho[1].TensorProduct(normalVector)) * dAd1Matrix;
			var mat2 = dAd1Matrix.Transpose() * (normalVector.TensorProduct(dRho[2])) * aMatrix + aMatrix.Transpose() * (dRho[1].TensorProduct(normalVector)) * dAd2Matrix;
			var mat3 = dAd2Matrix.Transpose() * (normalVector.TensorProduct(dRho[1])) * aMatrix + aMatrix.Transpose() * (dRho[2].TensorProduct(normalVector)) * dAd1Matrix;
			var mat4 = dAd2Matrix.Transpose() * (normalVector.TensorProduct(dRho[2])) * aMatrix + aMatrix.Transpose() * (dRho[2].TensorProduct(normalVector)) * dAd2Matrix;
			var Kr = mat1.Scale(scalar1) + mat2.Scale(scalar2) + mat3.Scale(scalar3) + mat4.Scale(scalar4);
			return Kr;
		}

		private Matrix CalculateTangentialStiffnessPartForSticking(Matrix aMatrix, Matrix da1Matrix, Matrix da2Matrix,
				Matrix mInv, double[] dRho1, double[] dRho2, double[] n, double[] T)
		{
			var surfaceVector1_x_surfaceVector1 = dRho1.TensorProduct(dRho1);
			var surfaceVector1_x_surfaceVector2 = dRho1.TensorProduct(dRho2);
			var surfaceVector2_x_surfaceVector1 = dRho2.TensorProduct(dRho1);
			var surfaceVector2_x_surfaceVector2 = dRho2.TensorProduct(dRho2);

			var TangentialStiffnessPart1 = (aMatrix.Transpose() * surfaceVector1_x_surfaceVector1 * aMatrix).Scale(PenaltyFactorTangential * mInv[0, 0]) +
				(aMatrix.Transpose() * surfaceVector1_x_surfaceVector2 * aMatrix).Scale(PenaltyFactorTangential * mInv[0, 1]) +
				(aMatrix.Transpose() * surfaceVector2_x_surfaceVector1 * aMatrix).Scale(PenaltyFactorTangential * mInv[1, 0]) +
				(aMatrix.Transpose() * surfaceVector2_x_surfaceVector2 * aMatrix).Scale(PenaltyFactorTangential * mInv[1, 1]);

			var TangentialStiffnessPart21 = (aMatrix.Transpose() * surfaceVector1_x_surfaceVector1 * da1Matrix).Scale(T[0] * mInv[0, 0] * mInv[0, 0]) +
				(da1Matrix.Transpose() * surfaceVector1_x_surfaceVector1 * aMatrix).Scale(T[0] * mInv[0, 0] * mInv[0, 0]);

			var TangentialStiffnessPart22 = (aMatrix.Transpose() * surfaceVector2_x_surfaceVector1 * da1Matrix).Scale(T[0] * mInv[0, 0] * mInv[0, 1]) +
				(da1Matrix.Transpose() * surfaceVector1_x_surfaceVector2 * aMatrix).Scale(T[0] * mInv[0, 0] * mInv[0, 1]);

			var TangentialStiffnessPart23 = (aMatrix.Transpose() * surfaceVector1_x_surfaceVector2 * da1Matrix).Scale(T[0] * mInv[0, 1] * mInv[0, 0]) +
				(da1Matrix.Transpose() * surfaceVector2_x_surfaceVector1 * aMatrix).Scale(T[0] * mInv[0, 1] * mInv[0, 0]);

			var TangentialStiffnessPart24 = (aMatrix.Transpose() * surfaceVector2_x_surfaceVector2 * da1Matrix).Scale(T[0] * mInv[0, 1] * mInv[0, 1]) +
				(da1Matrix.Transpose() * surfaceVector2_x_surfaceVector2 * aMatrix).Scale(T[0] * mInv[0, 1] * mInv[0, 1]);

			var TangentialStiffnessPart25 = (aMatrix.Transpose() * surfaceVector1_x_surfaceVector1 * da2Matrix).Scale(T[0] * mInv[0, 0] * mInv[1, 0]) +
				(da2Matrix.Transpose() * surfaceVector1_x_surfaceVector1 * aMatrix).Scale(T[0] * mInv[0, 0] * mInv[1, 0]);

			var TangentialStiffnessPart26 = (aMatrix.Transpose() * surfaceVector2_x_surfaceVector1 * da2Matrix).Scale(T[0] * mInv[0, 0] * mInv[1, 1]) +
				(da2Matrix.Transpose() * surfaceVector1_x_surfaceVector2 * aMatrix).Scale(T[0] * mInv[0, 0] * mInv[1, 1]);

			var TangentialStiffnessPart27 = (aMatrix.Transpose() * surfaceVector1_x_surfaceVector2 * da2Matrix).Scale(T[0] * mInv[0, 1] * mInv[1, 0]) +
				(da2Matrix.Transpose() * surfaceVector2_x_surfaceVector1 * aMatrix).Scale(T[0] * mInv[0, 1] * mInv[1, 0]);

			var TangentialStiffnessPart28 = (aMatrix.Transpose() * surfaceVector2_x_surfaceVector2 * da2Matrix).Scale(T[0] * mInv[0, 1] * mInv[1, 1]) +
				(da2Matrix.Transpose() * surfaceVector2_x_surfaceVector2 * aMatrix).Scale(T[0] * mInv[0, 1] * mInv[1, 1]);

			var TangentialStiffnessPart29 = (aMatrix.Transpose() * surfaceVector1_x_surfaceVector1 * da1Matrix).Scale(T[1] * mInv[1, 0] * mInv[0, 0]) +
				(da1Matrix.Transpose() * surfaceVector1_x_surfaceVector1 * aMatrix).Scale(T[1] * mInv[1, 0] * mInv[0, 0]);

			var TangentialStiffnessPart210 = (aMatrix.Transpose() * surfaceVector2_x_surfaceVector1 * da1Matrix).Scale(T[1] * mInv[1, 0] * mInv[0, 1]) +
				(da1Matrix.Transpose() * surfaceVector1_x_surfaceVector2 * aMatrix).Scale(T[1] * mInv[1, 0] * mInv[0, 1]);

			var TangentialStiffnessPart211 = (aMatrix.Transpose() * surfaceVector1_x_surfaceVector2 * da1Matrix).Scale(T[1] * mInv[1, 1] * mInv[0, 0]) +
				(da1Matrix.Transpose() * surfaceVector2_x_surfaceVector1 * aMatrix).Scale(T[1] * mInv[1, 1] * mInv[0, 0]);

			var TangentialStiffnessPart212 = (aMatrix.Transpose() * surfaceVector2_x_surfaceVector2 * da1Matrix).Scale(T[1] * mInv[1, 1] * mInv[0, 1]) +
				(da1Matrix.Transpose() * surfaceVector2_x_surfaceVector2 * aMatrix).Scale(T[1] * mInv[1, 1] * mInv[0, 1]);

			var TangentialStiffnessPart213 = (aMatrix.Transpose() * surfaceVector1_x_surfaceVector1 * da2Matrix).Scale(T[1] * mInv[1, 0] * mInv[1, 0]) +
				(da2Matrix.Transpose() * surfaceVector1_x_surfaceVector1 * aMatrix).Scale(T[1] * mInv[1, 0] * mInv[1, 0]);

			var TangentialStiffnessPart214 = (aMatrix.Transpose() * surfaceVector2_x_surfaceVector1 * da2Matrix).Scale(T[1] * mInv[1, 0] * mInv[1, 1]) +
				(da2Matrix.Transpose() * surfaceVector1_x_surfaceVector2 * aMatrix).Scale(T[1] * mInv[1, 0] * mInv[1, 1]);

			var TangentialStiffnessPart215 = (aMatrix.Transpose() * surfaceVector1_x_surfaceVector2 * da2Matrix).Scale(T[1] * mInv[1, 1] * mInv[1, 0]) +
				(da2Matrix.Transpose() * surfaceVector2_x_surfaceVector1 * aMatrix).Scale(T[1] * mInv[1, 1] * mInv[1, 0]);

			var TangentialStiffnessPart216 = (aMatrix.Transpose() * surfaceVector2_x_surfaceVector2 * da2Matrix).Scale(T[1] * mInv[1, 1] * mInv[1, 1]) +
				(da2Matrix.Transpose() * surfaceVector2_x_surfaceVector2 * aMatrix).Scale(T[1] * mInv[1, 1] * mInv[1, 1]);

			var TangentialStiffnessPart2 = TangentialStiffnessPart21 + TangentialStiffnessPart22 + TangentialStiffnessPart23 + TangentialStiffnessPart24 + TangentialStiffnessPart25 +
				TangentialStiffnessPart26 + TangentialStiffnessPart27 + TangentialStiffnessPart28 + TangentialStiffnessPart29 + TangentialStiffnessPart210 + TangentialStiffnessPart211 +
				TangentialStiffnessPart212 + TangentialStiffnessPart213 + TangentialStiffnessPart214 + TangentialStiffnessPart215 + TangentialStiffnessPart216;

			var TangentialStiffnessPart = (TangentialStiffnessPart1.Scale(-1d)) + TangentialStiffnessPart2.Scale(-1d);
			return TangentialStiffnessPart;
		}

		private Matrix CalculateTangentialStiffnessPartForSliding(Matrix aMatrix, Matrix da1Matrix, Matrix da2Matrix,
			Matrix mInv, double[] dRho1, double[] dRho2, double[] n,
			double[] T, double ksi3, double trialTangentialTractionNorm)
		{
			var surfaceVector1_x_surfaceVector1 = dRho1.TensorProduct(dRho1);
			var surfaceVector1_x_surfaceVector2 = dRho1.TensorProduct(dRho2);
			var surfaceVector2_x_surfaceVector1 = dRho2.TensorProduct(dRho1);
			var surfaceVector2_x_surfaceVector2 = dRho2.TensorProduct(dRho2);
			var surfaceVector1_x_n = dRho1.TensorProduct(n);
			var n_x_surfaceVector1 = n.TensorProduct(dRho1);
			var surfaceVector2_x_n = dRho2.TensorProduct(n);
			var n_x_surfaceVector2 = n.TensorProduct(dRho2);

			var TangentialStiffnessPart1 =
				((aMatrix.Transpose() * surfaceVector1_x_n * aMatrix).Scale(mInv[0, 0]) +
				(aMatrix.Transpose() * surfaceVector2_x_n * aMatrix).Scale(mInv[0, 1])).
				Scale(PenaltyFactorNormal * SlidingCoefficient * T[0] / trialTangentialTractionNorm)
				+
				((aMatrix.Transpose() * surfaceVector1_x_n * aMatrix).Scale(mInv[1, 0]) +
				(aMatrix.Transpose() * surfaceVector2_x_n * aMatrix).Scale(mInv[1, 1])).
				Scale(PenaltyFactorNormal * SlidingCoefficient * T[1] / trialTangentialTractionNorm);

			var TangentialStiffnessPart2 =
				((aMatrix.Transpose() * surfaceVector1_x_surfaceVector1 * aMatrix).Scale(mInv[0, 0]) +
				(aMatrix.Transpose() * surfaceVector1_x_surfaceVector2 * aMatrix).Scale(mInv[0, 1])).
				Scale(PenaltyFactorTangential * SlidingCoefficient * PenaltyFactorNormal * Math.Abs(ksi3) / trialTangentialTractionNorm)
				+
				((aMatrix.Transpose() * surfaceVector2_x_surfaceVector1 * aMatrix).Scale(mInv[1, 0]) +
				(aMatrix.Transpose() * surfaceVector2_x_surfaceVector2 * aMatrix).Scale(mInv[1, 1])).
				Scale(PenaltyFactorTangential * SlidingCoefficient * PenaltyFactorNormal * Math.Abs(ksi3) / trialTangentialTractionNorm);

			var matrix31 =
				aMatrix.Transpose() *
				(
				surfaceVector1_x_surfaceVector1.Scale(mInv[0, 0] * mInv[0, 0]) +
				surfaceVector1_x_surfaceVector2.Scale(mInv[0, 0] * mInv[0, 1]) +
				surfaceVector2_x_surfaceVector1.Scale(mInv[0, 1] * mInv[0, 0]) +
				surfaceVector2_x_surfaceVector2.Scale(mInv[0, 1] * mInv[0, 1])
				)
				* aMatrix;

			var matrix32 =
				aMatrix.Transpose() *
				(
				surfaceVector1_x_surfaceVector1.Scale(mInv[0, 0] * mInv[1, 0]) +
				surfaceVector1_x_surfaceVector2.Scale(mInv[0, 0] * mInv[1, 1]) +
				surfaceVector2_x_surfaceVector1.Scale(mInv[0, 1] * mInv[1, 0]) +
				surfaceVector2_x_surfaceVector2.Scale(mInv[0, 1] * mInv[1, 1])
				)
				* aMatrix;

			var matrix33 =
				aMatrix.Transpose() *
				(
				surfaceVector1_x_surfaceVector1.Scale(mInv[1, 0] * mInv[0, 0]) +
				surfaceVector1_x_surfaceVector2.Scale(mInv[1, 0] * mInv[0, 1]) +
				surfaceVector2_x_surfaceVector1.Scale(mInv[1, 1] * mInv[0, 0]) +
				surfaceVector2_x_surfaceVector2.Scale(mInv[1, 1] * mInv[0, 1])
				)
				* aMatrix;

			var matrix34 =
				aMatrix.Transpose() *
				(
				surfaceVector1_x_surfaceVector1.Scale(mInv[1, 0] * mInv[1, 0]) +
				surfaceVector1_x_surfaceVector2.Scale(mInv[1, 0] * mInv[1, 1]) +
				surfaceVector2_x_surfaceVector1.Scale(mInv[1, 1] * mInv[1, 0]) +
				surfaceVector2_x_surfaceVector2.Scale(mInv[1, 1] * mInv[1, 1])
				)
				* aMatrix;

			var TangentialStiffnessPart3 =
				matrix31.Scale
				(
				PenaltyFactorTangential * SlidingCoefficient * PenaltyFactorNormal * Math.Abs(ksi3) *
				T[0] * T[0] / Math.Pow(trialTangentialTractionNorm, 3)
				) +
				matrix32.Scale
				(
				PenaltyFactorTangential * SlidingCoefficient * PenaltyFactorNormal * Math.Abs(ksi3) *
				T[0] * T[1] / Math.Pow(trialTangentialTractionNorm, 3)
				) +
				matrix33.Scale
				(
				PenaltyFactorTangential * SlidingCoefficient * PenaltyFactorNormal * Math.Abs(ksi3) *
				T[1] * T[0] / Math.Pow(trialTangentialTractionNorm, 3)
				) +
				matrix34.Scale
				(
				PenaltyFactorTangential * SlidingCoefficient * PenaltyFactorNormal * Math.Abs(ksi3) *
				T[1] * T[1] / Math.Pow(trialTangentialTractionNorm, 3)
				);

			var matrix41 =
				aMatrix.Transpose() *
				(
				surfaceVector1_x_surfaceVector1.Scale(mInv[0, 0] * mInv[0, 0]) +
				surfaceVector2_x_surfaceVector1.Scale(mInv[0, 0] * mInv[0, 1]) +
				surfaceVector1_x_surfaceVector2.Scale(mInv[0, 1] * mInv[0, 0]) +
				surfaceVector2_x_surfaceVector2.Scale(mInv[0, 1] * mInv[0, 1])
				)
				* da1Matrix;

			var matrix42 =
				aMatrix.Transpose() *
				(
				surfaceVector1_x_surfaceVector1.Scale(mInv[0, 0] * mInv[1, 0]) +
				surfaceVector2_x_surfaceVector1.Scale(mInv[0, 0] * mInv[1, 1]) +
				surfaceVector1_x_surfaceVector2.Scale(mInv[0, 1] * mInv[1, 0]) +
				surfaceVector2_x_surfaceVector2.Scale(mInv[0, 1] * mInv[1, 1])
				)
				* da2Matrix;

			var matrix43 =
				da1Matrix.Transpose() *
				(
				surfaceVector1_x_surfaceVector1.Scale(mInv[0, 0] * mInv[0, 0]) +
				surfaceVector1_x_surfaceVector2.Scale(mInv[0, 0] * mInv[0, 1]) +
				surfaceVector2_x_surfaceVector1.Scale(mInv[0, 1] * mInv[0, 0]) +
				surfaceVector2_x_surfaceVector2.Scale(mInv[0, 1] * mInv[0, 1])
				)
				* aMatrix;

			var matrix44 =
				da2Matrix.Transpose() *
				(
				surfaceVector1_x_surfaceVector1.Scale(mInv[0, 0] * mInv[1, 0]) +
				surfaceVector1_x_surfaceVector2.Scale(mInv[0, 0] * mInv[1, 1]) +
				surfaceVector2_x_surfaceVector1.Scale(mInv[0, 1] * mInv[1, 0]) +
				surfaceVector2_x_surfaceVector2.Scale(mInv[0, 1] * mInv[1, 1])
				)
				* aMatrix;

			var matrix45 =
				aMatrix.Transpose() *
				(
				surfaceVector1_x_surfaceVector1.Scale(mInv[1, 0] * mInv[0, 0]) +
				surfaceVector2_x_surfaceVector1.Scale(mInv[1, 0] * mInv[0, 1]) +
				surfaceVector1_x_surfaceVector2.Scale(mInv[1, 1] * mInv[0, 0]) +
				surfaceVector2_x_surfaceVector2.Scale(mInv[1, 1] * mInv[0, 1])
				)
				* da1Matrix;

			var matrix46 =
				aMatrix.Transpose() *
				(
				surfaceVector1_x_surfaceVector1.Scale(mInv[1, 0] * mInv[1, 0]) +
				surfaceVector2_x_surfaceVector1.Scale(mInv[1, 0] * mInv[1, 1]) +
				surfaceVector1_x_surfaceVector2.Scale(mInv[1, 1] * mInv[1, 0]) +
				surfaceVector2_x_surfaceVector2.Scale(mInv[1, 1] * mInv[1, 1])
				)
				* da2Matrix;

			var matrix47 =
				da1Matrix.Transpose() *
				(
				surfaceVector1_x_surfaceVector1.Scale(mInv[1, 0] * mInv[0, 0]) +
				surfaceVector1_x_surfaceVector2.Scale(mInv[1, 0] * mInv[0, 1]) +
				surfaceVector2_x_surfaceVector1.Scale(mInv[1, 1] * mInv[0, 0]) +
				surfaceVector2_x_surfaceVector2.Scale(mInv[1, 1] * mInv[0, 1])
				)
				* aMatrix;

			var matrix48 =
				da2Matrix.Transpose() *
				(
				surfaceVector1_x_surfaceVector1.Scale(mInv[1, 0] * mInv[1, 0]) +
				surfaceVector1_x_surfaceVector2.Scale(mInv[1, 0] * mInv[1, 1]) +
				surfaceVector2_x_surfaceVector1.Scale(mInv[1, 1] * mInv[1, 0]) +
				surfaceVector2_x_surfaceVector2.Scale(mInv[1, 1] * mInv[1, 1])
				)
				* aMatrix;

			var TangentialStiffnessPart4 =
				(matrix41 + matrix42 + matrix43 + matrix44).Scale
				(
					SlidingCoefficient * PenaltyFactorNormal * Math.Abs(ksi3) * T[0] / trialTangentialTractionNorm
				)
				+
				(matrix45 + matrix46 + matrix47 + matrix48).Scale
				(
					SlidingCoefficient * PenaltyFactorNormal * Math.Abs(ksi3) * T[1] / trialTangentialTractionNorm
				);

			var TangentialStiffnessPart = TangentialStiffnessPart1 + TangentialStiffnessPart2 + TangentialStiffnessPart3.Scale(-1d) + TangentialStiffnessPart4;
			return TangentialStiffnessPart;
		}

		public IMatrix StiffnessMatrix()
		{
			var ksiVector = Project(new double[2]);
			if (Math.Abs(ksiVector[0]) <= 1.05 && Math.Abs(ksiVector[1]) <= 1.05)
			{
				var aMatrices = CalculatePositionMatrix(ksiVector[0], ksiVector[1]);
				var masterSurfaceVectors = SurfaceVectors(aMatrices.Item2, aMatrices.Item3);
				var m = MetricTensor(masterSurfaceVectors);
				var n = NormalVector(m, masterSurfaceVectors);
				var xUpdated = XUpdatedVector();
				var ksi3 = CalculatePenetration(n, aMatrices.Item1, xUpdated);
				if (ksi3 <= 0)
				{
					var Km = CalculateMainStiffnessPart(n, Matrix.CreateFromArray(aMatrices.Item1));
					var Kr = CalculateRotationalStiffnessPart(n, Matrix.CreateFromArray(aMatrices.Item1),
						Matrix.CreateFromArray(aMatrices.Item2), Matrix.CreateFromArray(aMatrices.Item3),
						masterSurfaceVectors, ksi3, m);
					var K = Km.Add(Kr);
					//
					var mInv = InverseMetricTensor(m);
					var xM = XUpdatedVector().Scale(-1d);
					var aMatrixOld = CalculatePositionMatrix(StickingPoint[1][0], StickingPoint[1][1]).Item1;
					var positionVector = new double[3];
					for (var k = 0; k < 4; k++)
					{
						positionVector[0] += xM[3 * k] * aMatrices.Item1[0, 3 * k];
						positionVector[1] += xM[3 * k + 1] * aMatrices.Item1[0, 3 * k];
						positionVector[2] += xM[3 * k + 2] * aMatrices.Item1[0, 3 * k];
					}
					var oldPositionVector = new double[3];
					var nodalDIsplacementsIncrements = CalculateNodalDisplacementsIncrements(xUpdated);
					for (var k = 0; k < 4; k++)
					{
						oldPositionVector[0] += -PreviousConvergedSolutionNodalCoordinates[3 * k] * aMatrixOld[0, 3 * k] - nodalDIsplacementsIncrements[3 * k] * aMatrixOld[0, 3 * k];
						oldPositionVector[1] += -PreviousConvergedSolutionNodalCoordinates[3 * k + 1] * aMatrixOld[0, 3 * k] - nodalDIsplacementsIncrements[3 * k + 1] * aMatrixOld[0, 3 * k];
						oldPositionVector[2] += -PreviousConvergedSolutionNodalCoordinates[3 * k + 2] * aMatrixOld[0, 3 * k] - nodalDIsplacementsIncrements[3 * k + 2] * aMatrixOld[0, 3 * k];
					}
					var deltaKsi = new double[]
					{
								(positionVector.Subtract(oldPositionVector)).DotProduct(masterSurfaceVectors[1]) * mInv[0,0] +
								(positionVector.Subtract(oldPositionVector)).DotProduct(masterSurfaceVectors[2]) * mInv[0,1],
								(positionVector.Subtract(oldPositionVector)).DotProduct(masterSurfaceVectors[1]) * mInv[1,0] +
								(positionVector.Subtract(oldPositionVector)).DotProduct(masterSurfaceVectors[2]) * mInv[1,1]
					};
					var rM1 = new double[]
					{
								CPPSurfaceBaseVector1[1][0],
								CPPSurfaceBaseVector1[1][1],
								CPPSurfaceBaseVector1[1][2]
					};
					var rM2 = new double[]
					{
								CPPSurfaceBaseVector2[1][0],
								CPPSurfaceBaseVector2[1][1],
								CPPSurfaceBaseVector2[1][2]
					};
					var mPrev = new double[,]
					{
								{ rM1.DotProduct(rM1), rM1.DotProduct(rM2) },
								{ rM2.DotProduct(rM1), rM2.DotProduct(rM2) }
					};
					var detmPrev = mPrev[0, 0] * mPrev[1, 1] - mPrev[0, 1] * mPrev[1, 0];
					var mPrevInv = new double[,]
					{
								{ mPrev[1,1]/detmPrev, - mPrev[0,1]/detmPrev },
								{ - mPrev[1,0]/detmPrev, mPrev[0,0]/detmPrev }
					};
					var trialTangentialTraction = new double[]
					{
								TangentialTraction[1][0] * mPrevInv[0, 0] * rM1.DotProduct(masterSurfaceVectors[1]) +
								TangentialTraction[1][0] * mPrevInv[0, 1] * rM2.DotProduct(masterSurfaceVectors[1]) +
								TangentialTraction[1][1] * mPrevInv[1, 0] * rM1.DotProduct(masterSurfaceVectors[1]) +
								TangentialTraction[1][1] * mPrevInv[1, 1] * rM2.DotProduct(masterSurfaceVectors[1]) -
								mInv[0, 0] * PenaltyFactorTangential * deltaKsi[0] - mInv[0, 1] * PenaltyFactorTangential * deltaKsi[1],
								TangentialTraction[1][0] * mPrevInv[0, 0] * rM1.DotProduct(masterSurfaceVectors[2]) +
								TangentialTraction[1][0] * mPrevInv[0, 1] * rM2.DotProduct(masterSurfaceVectors[2]) +
								TangentialTraction[1][1] * mPrevInv[1, 0] * rM1.DotProduct(masterSurfaceVectors[2]) +
								TangentialTraction[1][1] * mPrevInv[1, 1] * rM2.DotProduct(masterSurfaceVectors[2]) +
								-mInv[1, 0] * PenaltyFactorTangential * deltaKsi[0] - mInv[1, 1] * PenaltyFactorTangential * deltaKsi[1]
					};
					var tangentialTractionNorm = Math.Pow(mInv[0, 0] * trialTangentialTraction[0] * trialTangentialTraction[0] +
						mInv[0, 1] * trialTangentialTraction[0] * trialTangentialTraction[1] +
						mInv[1, 0] * trialTangentialTraction[1] * trialTangentialTraction[0] +
						mInv[1, 1] * trialTangentialTraction[1] * trialTangentialTraction[1],
						0.5);
					var phiTr = tangentialTractionNorm - StickingCoefficient * PenaltyFactorNormal * Math.Abs(ksi3);
					if (phiTr <= 0)
					{
						//stick
						var StifnessMatrixTangentialPart = CalculateTangentialStiffnessPartForSticking(Matrix.CreateFromArray(aMatrices.Item1),
							Matrix.CreateFromArray(aMatrices.Item2).Scale(-1d),
							Matrix.CreateFromArray(aMatrices.Item3).Scale(-1d),
							mInv, masterSurfaceVectors[1], masterSurfaceVectors[2], n, trialTangentialTraction);
						var stickStifnessMatrix = StifnessMatrixTangentialPart.Scale(-1d);
						K = K.Add(stickStifnessMatrix);
					}
					else
					{
						//slide
						var StifnessMatrixTangentialPart = CalculateTangentialStiffnessPartForSliding(Matrix.CreateFromArray(aMatrices.Item1), Matrix.CreateFromArray(aMatrices.Item2).Scale(-1d),
							Matrix.CreateFromArray(aMatrices.Item3).Scale(-1d), mInv, masterSurfaceVectors[1], masterSurfaceVectors[2],
							n, trialTangentialTraction, ksi3, tangentialTractionNorm);
							K = K.Add(StifnessMatrixTangentialPart);
					}
					return dofEnumerator.GetTransformedMatrix(K);
				}
				else
				{
					return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(new double[15, 15]));
				}
			}
			else
			{
				return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(new double[15, 15]));
			}
		}
		public IMatrix PhysicsMatrix() => StiffnessMatrix();
		public IMatrix MassMatrix()
		{
			double[,] massMatrix = new double[15, 15];
			return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(massMatrix));
		}

		public IMatrix DampingMatrix()
		{
			double[,] dampingMatrix = new double[15, 15];
			return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(dampingMatrix));
		}

		private double[] CreateInternalGlobalForcesVector()
		{
			var ksiVector = Project(new double[2]);
			if (Math.Abs(ksiVector[0]) <= 1.05 && Math.Abs(ksiVector[1]) <= 1.05)
			{
				var aMatrices = CalculatePositionMatrix(ksiVector[0], ksiVector[1]);
				var masterSurfaceVectors = SurfaceVectors(aMatrices.Item2, aMatrices.Item3);
				var m = MetricTensor(masterSurfaceVectors);
				var n = NormalVector(m, masterSurfaceVectors);
				var xUpdated = XUpdatedVector();
				var ksi3 = CalculatePenetration(n, aMatrices.Item1, xUpdated);
				if (ksi3 <= 0)
				{
					var AT = Matrix.CreateFromArray(aMatrices.Item1).Transpose();
					var internalGlobalForcesVector = AT.Multiply(n).Scale(PenaltyFactorNormal * ksi3);
					//
					var AT_rM1 = AT.Multiply(masterSurfaceVectors[1]);
					var AT_rM2 = AT.Multiply(masterSurfaceVectors[2]);
					var mInv = InverseMetricTensor(m);
					var xM = XUpdatedVector().Scale(-1d);
					var aMatrixOld = CalculatePositionMatrix(StickingPoint[1][0], StickingPoint[1][1]).Item1;
					var positionVector = new double[3];
					for (var k = 0; k < 4; k++)
					{
						positionVector[0] += xM[3 * k] * aMatrices.Item1[0, 3 * k];
						positionVector[1] += xM[3 * k + 1] * aMatrices.Item1[0, 3 * k];
						positionVector[2] += xM[3 * k + 2] * aMatrices.Item1[0, 3 * k];
					}
					var oldPositionVector = new double[3];
					var nodalDIsplacementsIncrements = CalculateNodalDisplacementsIncrements(xUpdated);
					for (var k = 0; k < 4; k++)
					{
						oldPositionVector[0] += -PreviousConvergedSolutionNodalCoordinates[3 * k] * aMatrixOld[0, 3 * k] - nodalDIsplacementsIncrements[3 * k] * aMatrixOld[0, 3 * k];
						oldPositionVector[1] += -PreviousConvergedSolutionNodalCoordinates[3 * k + 1] * aMatrixOld[0, 3 * k] - nodalDIsplacementsIncrements[3 * k + 1] * aMatrixOld[0, 3 * k];
						oldPositionVector[2] += -PreviousConvergedSolutionNodalCoordinates[3 * k + 2] * aMatrixOld[0, 3 * k] - nodalDIsplacementsIncrements[3 * k + 2] * aMatrixOld[0, 3 * k];
					}
					var deltaKsi = new double[]
					{
								(positionVector.Subtract(oldPositionVector)).DotProduct(masterSurfaceVectors[1]) * mInv[0,0] +
								(positionVector.Subtract(oldPositionVector)).DotProduct(masterSurfaceVectors[2]) * mInv[0,1],
								(positionVector.Subtract(oldPositionVector)).DotProduct(masterSurfaceVectors[1]) * mInv[1,0] +
								(positionVector.Subtract(oldPositionVector)).DotProduct(masterSurfaceVectors[2]) * mInv[1,1]
					};
					var rM1 = new double[]
					{
								CPPSurfaceBaseVector1[1][0],
								CPPSurfaceBaseVector1[1][1],
								CPPSurfaceBaseVector1[1][2]
					};
					var rM2 = new double[]
					{
								CPPSurfaceBaseVector2[1][0],
								CPPSurfaceBaseVector2[1][1],
								CPPSurfaceBaseVector2[1][2]
					};
					var mPrev = new double[,]
					{
								{ rM1.DotProduct(rM1), rM1.DotProduct(rM2) },
								{ rM2.DotProduct(rM1), rM2.DotProduct(rM2) }
					};
					var detmPrev = mPrev[0, 0] * mPrev[1, 1] - mPrev[0, 1] * mPrev[1, 0];
					var mPrevInv = new double[,]
					{
								{ mPrev[1,1]/detmPrev, - mPrev[0,1]/detmPrev },
								{ - mPrev[1,0]/detmPrev, mPrev[0,0]/detmPrev }
					};
					var trialTangentialTraction = new double[]
					{
								TangentialTraction[1][0] * mPrevInv[0, 0] * rM1.DotProduct(masterSurfaceVectors[1]) +
								TangentialTraction[1][0] * mPrevInv[0, 1] * rM2.DotProduct(masterSurfaceVectors[1]) +
								TangentialTraction[1][1] * mPrevInv[1, 0] * rM1.DotProduct(masterSurfaceVectors[1]) +
								TangentialTraction[1][1] * mPrevInv[1, 1] * rM2.DotProduct(masterSurfaceVectors[1]) -
								mInv[0, 0] * PenaltyFactorTangential * deltaKsi[0] - mInv[0, 1] * PenaltyFactorTangential * deltaKsi[1],
								TangentialTraction[1][0] * mPrevInv[0, 0] * rM1.DotProduct(masterSurfaceVectors[2]) +
								TangentialTraction[1][0] * mPrevInv[0, 1] * rM2.DotProduct(masterSurfaceVectors[2]) +
								TangentialTraction[1][1] * mPrevInv[1, 0] * rM1.DotProduct(masterSurfaceVectors[2]) +
								TangentialTraction[1][1] * mPrevInv[1, 1] * rM2.DotProduct(masterSurfaceVectors[2]) +
								-mInv[1, 0] * PenaltyFactorTangential * deltaKsi[0] - mInv[1, 1] * PenaltyFactorTangential * deltaKsi[1]
					};
					var tangentialTractionNorm = Math.Pow(mInv[0, 0] * trialTangentialTraction[0] * trialTangentialTraction[0] +
						mInv[0, 1] * trialTangentialTraction[0] * trialTangentialTraction[1] +
						mInv[1, 0] * trialTangentialTraction[1] * trialTangentialTraction[0] +
						mInv[1, 1] * trialTangentialTraction[1] * trialTangentialTraction[1],
						0.5);
					var phiTr = tangentialTractionNorm - StickingCoefficient * PenaltyFactorNormal * Math.Abs(ksi3);
					if (phiTr <= 0)
					{
						//stick
						var tangentialLocalForcesVector =
							(((AT_rM1.Scale(mInv[0, 0])).Add(AT_rM2.Scale(mInv[1, 0]))).Scale(-trialTangentialTraction[0])).
							Add(((AT_rM1.Scale(mInv[0, 1])).Add(AT_rM2.Scale(mInv[1, 1]))).Scale(-trialTangentialTraction[1]));
						internalGlobalForcesVector = internalGlobalForcesVector.Add(tangentialLocalForcesVector);
					}
					else
					{
						//slide
						var T = new double[]
						{
									(trialTangentialTraction[0] / tangentialTractionNorm) * SlidingCoefficient * PenaltyFactorNormal * Math.Abs(ksi3),
									(trialTangentialTraction[1] / tangentialTractionNorm) * SlidingCoefficient * PenaltyFactorNormal * Math.Abs(ksi3)
						};
						var tangentialLocalForcesVector =
							(((AT_rM1.Scale(mInv[0, 0])).Add(AT_rM2.Scale(mInv[1, 0]))).Scale(-T[0])).
							Add(((AT_rM1.Scale(mInv[0, 1])).Add(AT_rM2.Scale(mInv[1, 1]))).Scale(-T[1]));
						internalGlobalForcesVector = internalGlobalForcesVector.Add(tangentialLocalForcesVector);
					}
					return internalGlobalForcesVector;
				}
				else
				{
					var internalGlobalForcesVector = new double[15];
					return internalGlobalForcesVector;
				}
			}
			else
			{
				var internalGlobalForcesVector = new double[15];
				return internalGlobalForcesVector;
			}
		}

		public Tuple<double[], double[]> CalculateResponse(double[] local_Displacements)
		{
			//WARNING: 1) No strains are computed 2) localdDisplacements are not used.
			//double[] strains = null;
			////double[] stresses = null;

			//double[] forces = CalculateResponseIntegral(local_Displacements);
			//double[] stresses = Array.ConvertAll(forces, x => x / ContactArea);
			if (DisplacementVector == null || DisplacementVector.Length != local_Displacements.Length)
			{
				DisplacementVector = new double[local_Displacements.Length];
			}

			Array.Copy(local_Displacements, DisplacementVector, local_Displacements.Length);

			//return new Tuple<double[], double[]>(strains, stresses);
			return new Tuple<double[], double[]>(null, null);
		}

		public double[] CalculateResponseIntegralForLogging(double[] localDisplacements)
		{
			return CalculateResponseIntegral();
		}

		public double[] CalculateResponseIntegral() => CreateInternalGlobalForcesVector();

		public void SaveConstitutiveLawState(IHaveState externalState) => UpdateTangentialProperties();

		#endregion

		#region IFiniteElement Members


		public bool ConstitutiveLawModified => false;

		public void ResetConstitutiveLawModified() { }

		#endregion

		#region IFiniteElement Members

		public void ClearConstitutiveLawState() { }

		public void ClearConstitutiveLawStresses() => throw new NotImplementedException();

		public IEnumerable<double[]> InterpolateElementModelQuantities(IEnumerable<IElementModelQuantity<IStructuralDofType>> quantities) => throw new NotImplementedException();
		public IEnumerable<double[]> IntegrateElementModelQuantities(IEnumerable<IElementModelQuantity<IStructuralDofType>> quantities) => throw new NotImplementedException();
		#endregion
	}
}
