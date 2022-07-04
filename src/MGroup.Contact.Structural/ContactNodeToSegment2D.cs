using System;
using System.Collections.Generic;
using MGroup.Constitutive.Structural;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.FEM.Structural.Line
{
	public class ContactNodeToSegment2D : IStructuralElementType
	{
		private static readonly IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY };
		private static readonly IDofType[][] dofs = new IDofType[][] { nodalDOFTypes, nodalDOFTypes, nodalDOFTypes };
		private readonly double penaltyFactor;
		private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
		private double[] DisplacementVector { get; set; }
		private double ContactArea { get; }

		public ContactNodeToSegment2D(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplier,
			double contactArea)
		{
			this.penaltyFactor = penaltyFactorMultiplier * youngModulus;
			this.DisplacementVector = new double[6];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
		}
		public ContactNodeToSegment2D(IReadOnlyList<INode> nodes, double penaltyFactor, double contactArea)
		{
			this.penaltyFactor = penaltyFactor;
			this.DisplacementVector = new double[6];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
		}
		public ContactNodeToSegment2D(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplier,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, youngModulus, penaltyFactorMultiplier, 1d)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public ContactNodeToSegment2D(IReadOnlyList<INode> nodes, double penaltyFactor, IElementDofEnumerator dofEnumerator)
			: this(nodes, penaltyFactor, 1.0)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public CellType CellType { get; } = CellType.Line2;

		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		public double ClosestPointProjection()
		{
			var Xm1 = Nodes[0].X + DisplacementVector[0];
			var Ym1 = Nodes[0].Y + DisplacementVector[1];
			var Xm2 = Nodes[1].X + DisplacementVector[2];
			var Ym2 = Nodes[1].Y + DisplacementVector[3];
			var Xs = Nodes[2].X + DisplacementVector[4];
			var Ys = Nodes[2].Y + DisplacementVector[5];

			var ksi1 = (2.0 * (Xs * (Xm2 - Xm1) + Ys * (Ym2 - Ym1)) - Math.Pow(Xm2, 2) - Math.Pow(Ym2, 2) + Math.Pow(Xm1, 2) + Math.Pow(Ym1, 2)) / (Math.Pow(Xm2 - Xm1, 2) + Math.Pow(Ym2 - Ym1, 2));
			return ksi1;
		}

		private double CalculateNormalGap(double[,] aMatrix, double[] n)
		{
			var AT = Matrix.CreateFromArray(aMatrix).Transpose();
			var AT_n = Vector.CreateFromArray(AT.Multiply(n));
			var xupd = Vector.CreateFromArray( new double[]
			{
				Nodes[0].X + DisplacementVector[0],
				Nodes[0].Y + DisplacementVector[1],
				Nodes[1].X + DisplacementVector[2],
				Nodes[1].Y + DisplacementVector[3],
				Nodes[2].X + DisplacementVector[4],
				Nodes[2].Y + DisplacementVector[5]
			});
			var normalGap = xupd.DotProduct(AT_n);
			return normalGap;
		}

		private Tuple<double[], double, double[]> MasterSegmentGeoParameters(double[,] daMatrix)
		{
			var Xm1 = Nodes[0].X + DisplacementVector[0];
			var Ym1 = Nodes[0].Y + DisplacementVector[1];
			var Xm2 = Nodes[1].X + DisplacementVector[2];
			var Ym2 = Nodes[1].Y + DisplacementVector[3];
			var Xs = Nodes[2].X + DisplacementVector[4];
			var Ys = Nodes[2].Y + DisplacementVector[5];

			var xupd = new double[] { -Xm1, -Ym1, -Xm2, -Ym2, -Xs, -Ys };
			var surfaceVector = Matrix.CreateFromArray(daMatrix).Multiply(xupd);
			var detm = surfaceVector.DotProduct(surfaceVector);
			var m11 = 1.0 / detm;
			var vector = new double[] { Ym2 - Ym1, Xm1 - Xm2 };
			var scalarCoef = -1.0 / (2.0 * Math.Sqrt(detm));
			var normalUnitVec = vector.Scale(scalarCoef);

			return new Tuple<double[], double, double[]>(surfaceVector, m11, normalUnitVec);
		}

		private Tuple<double[,], double[,]> CalculatePositionMatrix(double ksi1)
		{
			var N1 = 1.0 / 2.0 * (1.0 - ksi1);
			var N2 = 1.0 / 2.0 * (1.0 + ksi1);
			var dN1 = -1.0 / 2.0;
			var dN2 = 1.0 / 2.0;
			var aMatrix = new double[,]
				{
					{ -N1 ,0.0 ,-N2 ,0.0 ,1.0 ,0.0 },
					{0.0, -N1 , 0.0 ,-N2, 0.0, 1.0 }
				};

			var daMatrix = new double[,]
				{
					{ -dN1 ,0.0 ,-dN2 ,0.0 ,0.0 ,0.0 },
					{0.0, -dN1 , 0.0 ,-dN2, 0.0, 0.0 }
				};
			return new Tuple<double[,], double[,]>(aMatrix, daMatrix);
		}

		#region IElementType Members

		public int ID { get; set; }
		public int SubdomainID { get; set; }
		public IReadOnlyList<INode> Nodes { get; }

		public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofs;

		public IReadOnlyList<INode> GetNodesForMatrixAssembly() => Nodes;

		private Matrix CalculateMainStiffnessPart(double[] n, Matrix aMatrix)
		{
			var nxn = n.TensorProduct(n);
			var nxn_A = nxn.MultiplyRight(aMatrix);
			var AT_nxn_A = aMatrix.Transpose().MultiplyRight(nxn_A);
			var mainStiffnessMatrix = AT_nxn_A.Scale(penaltyFactor);
			return mainStiffnessMatrix;
		}

		private Matrix CalculateRotationalStiffnessPart(Matrix A, Matrix dA, double[] n, double penetration, double m11, double[] surfaceVector)
		{
			var coef = penaltyFactor * penetration * m11;
			var n_x_dRho = n.TensorProduct(surfaceVector);
			var dRho_x_n = surfaceVector.TensorProduct(n);
			var firstTerm = dA.Transpose().MultiplyRight(n_x_dRho.MultiplyRight(A));
			var secondTerm = A.Transpose().MultiplyRight(dRho_x_n.MultiplyRight(dA));
			var rotationalPart = (firstTerm + secondTerm).Scale(coef);
			return rotationalPart;
		}

		public IMatrix StiffnessMatrix()
		{
			var ksi1 = ClosestPointProjection();
			if (Math.Abs(ksi1) <= 1.05)
			{
				var positionMatrices = CalculatePositionMatrix(ksi1);
				var A = positionMatrices.Item1;
				var dA = positionMatrices.Item2;
				var masterSegmentGeometry = MasterSegmentGeoParameters(dA);
				var surfaceVector = masterSegmentGeometry.Item1;
				var metricTensorContravariant = masterSegmentGeometry.Item2;
				var n = masterSegmentGeometry.Item3;
				var penatration = CalculateNormalGap(A,n);
				if (penatration <= 0)
				{
					var mainPart = CalculateMainStiffnessPart(n, Matrix.CreateFromArray(A));
					var rotationalPart = CalculateRotationalStiffnessPart(Matrix.CreateFromArray(A), Matrix.CreateFromArray(dA),
										n, penatration, metricTensorContravariant, surfaceVector);
					var globalStiffnessMatrix = mainPart.Add(rotationalPart);
					return dofEnumerator.GetTransformedMatrix(globalStiffnessMatrix);
				}
				else
				{
					var globalStifnessMatrix = new double[6, 6];
					return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(globalStifnessMatrix));
				}

			}
			else
			{
				var globalStifnessMatrix = new double[6, 6];
				return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(globalStifnessMatrix));
			}
		}

		public IMatrix PhysicsMatrix() => StiffnessMatrix();

		public IMatrix MassMatrix()
		{
			double[,] massMatrix = new double[6, 6];
			return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(massMatrix));
		}

		public IMatrix DampingMatrix()
		{
			double[,] dampingMatrix = new double[6, 6];
			return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(dampingMatrix));
		}

		private double[] CreateInternalGlobalForcesVector()
		{
			var ksi1 = ClosestPointProjection();
			if (Math.Abs(ksi1) <= 1.05)
			{
				var positionMatrices = CalculatePositionMatrix(ksi1);
				var A = positionMatrices.Item1;
				var dA = positionMatrices.Item2;
				var masterSegmentGeometry = MasterSegmentGeoParameters(dA);
				var n = masterSegmentGeometry.Item3;
				var penetration = CalculateNormalGap(A, n);
				if (penetration <= 0)
				{
					var AT = Matrix.CreateFromArray(A).Transpose();
					var AT_n = AT.Multiply(n);
					var internalGlobalForcesVector = AT_n.Scale(penaltyFactor * penetration);
					return internalGlobalForcesVector;
				}
				else
				{
					var internalGlobalForcesVector = new double[6];
					return internalGlobalForcesVector;
				}
			}
			else
			{
				var internalGlobalForcesVector = new double[6];
				return internalGlobalForcesVector;
			}
		}

		public Tuple<double[], double[]> CalculateResponse(double[] local_Displacements)
		{
			// WARNING: 1) No strains are computed 2) localdDisplacements are not used.
			double[] strains = null;
			double[] forces = CreateInternalGlobalForcesVector();
			double[] stresses = Array.ConvertAll(forces, x => x / ContactArea);
			if (DisplacementVector == null || DisplacementVector.Length != local_Displacements.Length)
			{
				DisplacementVector = new double[local_Displacements.Length];
			}

			Array.Copy(local_Displacements, DisplacementVector, local_Displacements.Length);

			return new Tuple<double[], double[]>(strains, stresses);
		}

		public double[] CalculateResponseIntegralForLogging(double[] localDisplacements) => CalculateResponseIntegral();

		public double[] CalculateResponseIntegral() => CreateInternalGlobalForcesVector();

		public void SaveConstitutiveLawState() { }

		#endregion

		#region IFiniteElement Members


		public bool ConstitutiveLawModified => false;

		public void ResetConstitutiveLawModified() { }

		#endregion

		#region IFiniteElement Members

		public void ClearConstitutiveLawState() { }

		public void ClearConstitutiveLawStresses() => throw new NotImplementedException();

		#endregion
	}
}
