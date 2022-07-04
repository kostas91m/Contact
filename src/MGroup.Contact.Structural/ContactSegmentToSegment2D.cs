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
	public class ContactSegmentToSegment2D : IStructuralElementType
	{
		private readonly IDofType[][] dofs;
		private readonly double penaltyFactor;
		private readonly int MasterSegmentOrder;
		private readonly int SlaveSegmentOrder;
		private readonly int IntegrationPoints;

		private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
		private double[] DisplacementVector { get; set; }
		private double ContactArea { get; }

		public ContactSegmentToSegment2D(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplier, double contactArea)
		{
			this.penaltyFactor = penaltyFactorMultiplier * youngModulus;
			this.DisplacementVector = new double[2 * nodes.Count];
			this.ContactArea = contactArea;
			if (nodes.Count != 4)
			{
				throw new ArgumentException("This Constructor can be used only for linear segments");
			}
			else
			{
				this.Nodes = nodes;
				dofs = new IDofType[nodes.Count][];
				for (var i = 0; i < nodes.Count; i++)
				{
					dofs[i] = new IDofType[]
					{
						StructuralDof.TranslationX, StructuralDof.TranslationY
					};
				}
				this.MasterSegmentOrder = 1;
				this.SlaveSegmentOrder = 1;
				this.IntegrationPoints = 2;
			}
		}
		public ContactSegmentToSegment2D(IReadOnlyList<INode> nodes, double penaltyFactor, double contactArea)
		{
			this.penaltyFactor = penaltyFactor;
			this.DisplacementVector = new double[2 * nodes.Count];
			this.ContactArea = contactArea;
			if (nodes.Count != 4)
			{
				throw new ArgumentException("This Constructor can be used only for linear segments");
			}
			else
			{
				this.Nodes = nodes;
				dofs = new IDofType[nodes.Count][];
				for (var i = 0; i < nodes.Count; i++)
				{
					dofs[i] = new IDofType[]
					{
						StructuralDof.TranslationX, StructuralDof.TranslationY
					};
				}
				this.MasterSegmentOrder = 1;
				this.SlaveSegmentOrder = 1;
				this.IntegrationPoints = 2;
			}
		}
		public ContactSegmentToSegment2D(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplier, double contactArea,
			int masterSegmentOrder, int slaveSegmentOrder)
		{
			if((masterSegmentOrder!= 1 && masterSegmentOrder != 2) || (slaveSegmentOrder != 1 && slaveSegmentOrder != 2))
			{
				throw new ArgumentException("Only linear & quadratic segments can be defined");
			}
			else
			{
				this.penaltyFactor = penaltyFactorMultiplier * youngModulus;
				this.DisplacementVector = new double[2 * nodes.Count];
				this.ContactArea = contactArea;
				if (nodes.Count != slaveSegmentOrder + masterSegmentOrder + 2)
				{
					throw new ArgumentException("Inconsistent input regarding the nodes & the order of the segments");
				}
				else
				{
					this.Nodes = nodes;
					dofs = new IDofType[nodes.Count][];
					for (var i = 0; i < nodes.Count; i++)
					{
						dofs[i] = new IDofType[]
						{
							StructuralDof.TranslationX, StructuralDof.TranslationY
						};
					}
					this.MasterSegmentOrder = masterSegmentOrder;
					this.SlaveSegmentOrder = slaveSegmentOrder;
					if (slaveSegmentOrder == 1)
					{
						this.IntegrationPoints = 2;
					}
					else
					{
						this.IntegrationPoints = 3;
					}
				}
			}			
		}
		public ContactSegmentToSegment2D(IReadOnlyList<INode> nodes, double penaltyFactor, double contactArea, int masterSegmentOrder, int slaveSegmentOrder)
		{
			if ((masterSegmentOrder != 1 && masterSegmentOrder != 2) || (slaveSegmentOrder != 1 && slaveSegmentOrder != 2))
			{
				throw new ArgumentException("Only linear & quadratic segments can be defined");
			}
			else
			{
				this.penaltyFactor = penaltyFactor;
				this.DisplacementVector = new double[2 * nodes.Count];
				this.ContactArea = contactArea;
				if (nodes.Count != slaveSegmentOrder + masterSegmentOrder + 2)
				{
					throw new ArgumentException("Inconsistent input regarding the nodes & the order of the segments");
				}
				else
				{
					this.Nodes = nodes;
					dofs = new IDofType[nodes.Count][];
					for (var i = 0; i < nodes.Count; i++)
					{
						dofs[i] = new IDofType[]
						{
							StructuralDof.TranslationX, StructuralDof.TranslationY
						};
					}
					this.MasterSegmentOrder = masterSegmentOrder;
					this.SlaveSegmentOrder = slaveSegmentOrder;
					if(slaveSegmentOrder == 1)
					{
						this.IntegrationPoints = 2;
					}
					else
					{
						this.IntegrationPoints = 3;
					}
				}
			}
		}
		public ContactSegmentToSegment2D(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplier, double contactArea,
	int masterSegmentOrder, int slaveSegmentOrder, int integrationPoints)
		{
			if ((masterSegmentOrder != 1 && masterSegmentOrder != 2) || (slaveSegmentOrder != 1 && slaveSegmentOrder != 2))
			{
				throw new ArgumentException("Only linear & quadratic segments can be defined");
			}
			else
			{
				this.penaltyFactor = penaltyFactorMultiplier * youngModulus;
				this.DisplacementVector = new double[2 * nodes.Count];
				this.ContactArea = contactArea;
				if (nodes.Count != slaveSegmentOrder + masterSegmentOrder + 2)
				{
					throw new ArgumentException("Inconsistent input regarding the nodes & the order of the segments");
				}
				else
				{
					this.Nodes = nodes;
					dofs = new IDofType[nodes.Count][];
					for (var i = 0; i < nodes.Count; i++)
					{
						dofs[i] = new IDofType[]
						{
							StructuralDof.TranslationX, StructuralDof.TranslationY
						};
					}
					this.MasterSegmentOrder = masterSegmentOrder;
					this.SlaveSegmentOrder = slaveSegmentOrder;
					if(integrationPoints < 2 || integrationPoints > 6)
					{
						throw new ArgumentException("Between [2,6] Gauss points can be defined");
					}
					this.IntegrationPoints = integrationPoints;
				}
			}
		}
		public ContactSegmentToSegment2D(IReadOnlyList<INode> nodes, double penaltyFactor, double contactArea, int masterSegmentOrder, int slaveSegmentOrder, int integrationPoints)
		{
			if ((masterSegmentOrder != 1 && masterSegmentOrder != 2) || (slaveSegmentOrder != 1 && slaveSegmentOrder != 2))
			{
				throw new ArgumentException("Only linear & quadratic segments can be defined");
			}
			else
			{
				this.penaltyFactor = penaltyFactor;
				this.DisplacementVector = new double[2 * nodes.Count];
				this.ContactArea = contactArea;
				if (nodes.Count != slaveSegmentOrder + masterSegmentOrder + 2)
				{
					throw new ArgumentException("Inconsistent input regarding the nodes & the order of the segments");
				}
				else
				{
					this.Nodes = nodes;
					dofs = new IDofType[nodes.Count][];
					for (var i = 0; i < nodes.Count; i++)
					{
						dofs[i] = new IDofType[]
						{
							StructuralDof.TranslationX, StructuralDof.TranslationY
						};
					}
					this.MasterSegmentOrder = masterSegmentOrder;
					this.SlaveSegmentOrder = slaveSegmentOrder;
					if (integrationPoints < 2 || integrationPoints > 6)
					{
						throw new ArgumentException("Between [2,6] Gauss points can be defined");
					}
					this.IntegrationPoints = integrationPoints;
				}
			}
		}
		public ContactSegmentToSegment2D(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplier,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, youngModulus, penaltyFactorMultiplier, 1d)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public ContactSegmentToSegment2D(IReadOnlyList<INode> nodes, double penaltyFactor, IElementDofEnumerator dofEnumerator)
			: this(nodes, penaltyFactor, 1.0)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public ContactSegmentToSegment2D(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplier, int masterSegmentOrder, int slaveSegmentOrder,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, youngModulus, penaltyFactorMultiplier, 1d, masterSegmentOrder, slaveSegmentOrder)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public ContactSegmentToSegment2D(IReadOnlyList<INode> nodes, double penaltyFactor, int masterSegmentOrder, int slaveSegmentOrder,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, penaltyFactor, 1.0, masterSegmentOrder, slaveSegmentOrder)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public ContactSegmentToSegment2D(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplier, int masterSegmentOrder, int slaveSegmentOrder, int integrationPoints,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, youngModulus, penaltyFactorMultiplier, 1d, masterSegmentOrder, slaveSegmentOrder, integrationPoints)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public ContactSegmentToSegment2D(IReadOnlyList<INode> nodes, double penaltyFactor, int masterSegmentOrder, int slaveSegmentOrder, int integrationPoints,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, penaltyFactor, 1.0, masterSegmentOrder, slaveSegmentOrder, integrationPoints)
		{
			this.dofEnumerator = dofEnumerator;
		}

		public CellType CellType { get; } = CellType.Line2;

		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}

		private double[] NodalXUpdated()
		{
			var x = new double[2 * Nodes.Count];
			for (var i = 0; i < Nodes.Count; i++)
			{
				x[i * 2] =	Nodes[i].X + DisplacementVector[i * 2];
				x[i * 2 + 1] = Nodes[i].Y + DisplacementVector[i * 2 + 1];
			}
			return x;
		}

		private Tuple<double[,], double[,], double[,], double[,], double[,]> CalculatePositionMatrix(double ksi1, double ksi2)
		{
			if (SlaveSegmentOrder == 1 && MasterSegmentOrder == 1)
			{
				var N1 = 1.0 / 2.0 * (1.0 - ksi1);
				var N2 = 1.0 / 2.0 * (1.0 + ksi1);
				var N3 = 1.0 / 2.0 * (1.0 - ksi2);
				var N4 = 1.0 / 2.0 * (1.0 + ksi2);
				var dN11 = -1.0 / 2.0;
				var dN21 = 1.0 / 2.0;
				var dN32 = -1.0 / 2.0;
				var dN42 = 1.0 / 2.0;
				var aMatrix = new double[,]
				{
					{ -N1, 0.0, -N2, 0.0, N3, 0.0, N4, 0.0 },
					{ 0.0, -N1, 0.0, -N2, 0.0, N3, 0.0, N4 }
				};
				var da1Matrix = new double[,]
				{
					{ -dN11, 0.0, -dN21, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN11, 0.0, -dN21, 0.0, 0.0, 0.0, 0.0 }
				};
				var da11Matrix = new double[,]
				{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
				};
				var da2Matrix = new double[,]
				{
					{ 0.0, 0.0, 0.0, 0.0, dN32, 0.0, dN42, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, dN32, 0.0, dN42 }
				};
				var da22Matrix = new double[,]
				{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
				};
				return new Tuple<double[,], double[,], double[,], double[,], double[,]>(aMatrix, da1Matrix, da11Matrix, da2Matrix, da22Matrix);
			}
			else if (MasterSegmentOrder == 1 && SlaveSegmentOrder == 2)
			{
				var N1 = 1.0 / 2.0 * (1.0 - ksi1);
				var N2 = 1.0 / 2.0 * (1.0 + ksi1);
				var N3 = 1.0 / 2.0 * (Math.Pow(ksi2, 2.0) - ksi2);
				var N4 = 1.0 - Math.Pow(ksi2, 2.0);
				var N5 = 1.0 / 2.0 * (Math.Pow(ksi2, 2.0) + ksi2);
				var dN11 = -1.0 / 2.0;
				var dN21 = 1.0 / 2.0;
				var dN32 = ksi2 - 0.5;
				var dN42 = -2.0 * ksi2;
				var dN52 = ksi2 + 0.5;
				var dN322 = 1.0;
				var dN422 = -2.0;
				var dN522 = 1.0;
				var aMatrix = new double[,]
				{
					{ -N1, 0.0, -N2, 0.0, N3, 0.0 , N4, 0.0, N5, 0.0 },
					{ 0.0, -N1, 0.0, -N2, 0.0, N3, 0.0, N4, 0.0, N5 }
				};
				var da1Matrix = new double[,]
				{
					{ -dN11, 0.0, -dN21, 0.0, 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN11, 0.0, -dN21, 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0 }
				};
				var da11Matrix = new double[,]
				{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0 }
				};
				var da2Matrix = new double[,]
				{
					{ 0.0, 0.0, 0.0, 0.0, dN32, 0.0 , dN42, 0.0, dN52, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, dN32 , 0.0, dN42, 0.0, dN52 }
				};
				var da22Matrix = new double[,]
				{
					{ 0.0, 0.0, 0.0, 0.0, dN322, 0.0 , dN422, 0.0, dN522, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, dN322 , 0.0, dN422, 0.0, dN522 }
				};
				return new Tuple<double[,], double[,], double[,], double[,], double[,]>(aMatrix, da1Matrix, da11Matrix, da2Matrix, da22Matrix);
			}
			else if (MasterSegmentOrder == 2 && SlaveSegmentOrder == 1)
			{

				var N1 = 1.0 / 2.0 * (Math.Pow(ksi1, 2.0) - ksi1);
				var N2 = 1.0 - Math.Pow(ksi1, 2.0);
				var N3 = 1.0 / 2.0 * (Math.Pow(ksi1, 2.0) + ksi1);
				var N4 = 1.0 / 2.0 * (1.0 - ksi2);
				var N5 = 1.0 / 2.0 * (1.0 + ksi2);
				var dN11 = ksi1 - 0.5;
				var dN21 = -2.0 * ksi1;
				var dN31 = ksi1 + 0.5;
				var dN42 = -1.0 / 2.0;
				var dN52 = 1.0 / 2.0;
				var dN111 = 1.0;
				var dN211 = -2.0;
				var dN311 = 1.0;
				var aMatrix = new double[,]
				{
					{ -N1 ,0.0 ,-N2 ,0.0 ,-N3 , 0.0 , N4, 0.0, N5, 0.0 },
					{0.0, -N1 , 0.0 ,-N2, 0.0, -N3, 0.0, N4, 0.0, N5 }
				};
				var da1Matrix = new double[,]
				{
					{ -dN11 ,0.0 ,-dN21 ,0.0 ,-dN31 , 0.0 , 0.0, 0.0, 0.0, 0.0 },
					{0.0, -dN11 , 0.0 ,-dN21, 0.0, -dN31, 0.0, 0.0, 0.0, 0.0 }
				};
				var da11Matrix = new double[,]
				{
					{ -dN111 ,0.0 ,-dN211 ,0.0 ,-dN311 , 0.0 , 0.0, 0.0, 0.0, 0.0 },
					{0.0, -dN111 , 0.0 ,-dN211, 0.0, -dN311, 0.0, 0.0, 0.0, 0.0 }
				};
				var da2Matrix = new double[,]
				{
					{ 0.0, 0.0, 0.0 ,0.0 ,0.0 , 0.0 , dN42, 0.0, dN52, 0.0 },
					{0.0, 0.0 , 0.0 ,0.0, 0.0, 0.0, 0.0, dN42, 0.0, dN52 }
				};
				var da22Matrix = new double[,]
				{
					{ 0.0, 0.0 ,0.0 ,0.0 ,0.0 , 0.0 , 0.0, 0.0, 0.0, 0.0 },
					{0.0, 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
				};
				return new Tuple<double[,], double[,], double[,], double[,], double[,]>(aMatrix, da1Matrix, da11Matrix, da2Matrix, da22Matrix);
			}
			else
			{
				var N1 = 1.0 / 2.0 * (Math.Pow(ksi1, 2.0) - ksi1);
				var N2 = 1.0 - Math.Pow(ksi1, 2.0);
				var N3 = 1.0 / 2.0 * (Math.Pow(ksi1, 2.0) + ksi1);
				var N4 = 1.0 / 2.0 * (Math.Pow(ksi2, 2.0) - ksi2);
				var N5 = 1.0 - Math.Pow(ksi2, 2.0);
				var N6 = 1.0 / 2.0 * (Math.Pow(ksi2, 2.0) + ksi2);
				var dN11 = ksi1 - 0.5;
				var dN21 = -2.0 * ksi1;
				var dN31 = ksi1 + 0.5;
				var dN42 = ksi2 - 0.5;
				var dN52 = -2.0 * ksi2;
				var dN62 = ksi2 + 0.5;
				var dN111 = 1.0;
				var dN211 = -2.0;
				var dN311 = 1.0;
				var dN422 = 1.0;
				var dN522 = -2.0;
				var dN622 = 1.0;
				var aMatrix = new double[,]
				{
					{ -N1 ,0.0 ,-N2 ,0.0 ,-N3 , 0.0 , N4, 0.0, N5, 0.0, N6, 0.0 },
					{0.0, -N1 , 0.0 ,-N2, 0.0, -N3, 0.0, N4, 0.0, N5, 0.0, N6 }
				};
				var da1Matrix = new double[,]
				{
					{ -dN11 ,0.0 ,-dN21 ,0.0 ,-dN31, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{0.0, -dN11 , 0.0 ,-dN21, 0.0, -dN31, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
				};
				var da11Matrix = new double[,]
				{
					{ -dN111 ,0.0 ,-dN211 ,0.0 ,-dN311, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{0.0, -dN111 , 0.0 ,-dN211, 0.0, -dN311, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
				};
				var da2Matrix = new double[,]
				{
					{ 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, dN42, 0.0, dN52, 0.0, dN62, 0.0 },
					{ 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, dN42, 0.0, dN52, 0.0, dN62 }
				};
				var da22Matrix = new double[,]
				{
					{ 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, dN422, 0.0, dN522, 0.0, dN622, 0.0 },
					{ 0.0, 0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, dN422, 0.0, dN522, 0.0, dN622 }
				};
				return new Tuple<double[,], double[,], double[,], double[,], double[,]>(aMatrix, da1Matrix, da11Matrix, da2Matrix, da22Matrix);
			}
		}

		private Tuple<double[], double, double[], double[], double> MasterSegmentGeometry(double[,] daMatrix, double[,] da2Matrix)
		{
			var xupd = NodalXUpdated().Scale(-1.0);
			var surfaceVector = Matrix.CreateFromArray(daMatrix).Multiply(xupd);
			var surfaceVectorDerivative = Matrix.CreateFromArray(da2Matrix).Multiply(xupd);
			var detm = surfaceVector.DotProduct(surfaceVector);
			var m11 = 1.0 / detm;
			var vector = new double[2];
			var scalarCoef = new double();
			var scalarCoef2 = Math.Pow(m11, 2.0);
			var tangentVector = new double[2];
			var curvatureTensor = new double();

			if (MasterSegmentOrder == 1)
			{
				var Xm1 = Nodes[0].X + DisplacementVector[0];
				var Ym1 = Nodes[0].Y + DisplacementVector[1];
				var Xm2 = Nodes[1].X + DisplacementVector[2];
				var Ym2 = Nodes[1].Y + DisplacementVector[3];
				vector[0] = Ym2 - Ym1;
				vector[1] = Xm1 - Xm2;
				scalarCoef = -1.0 / (2.0 * Math.Sqrt(detm));
			}
			else
			{
				vector[0] = -surfaceVector[1];
				vector[1] = surfaceVector[0];
				scalarCoef = 1.0 / (Math.Sqrt(detm));
				tangentVector = surfaceVector.Scale(scalarCoef);
			}
			var normalUnitVec = vector.Scale(scalarCoef);
			curvatureTensor = scalarCoef2 * surfaceVectorDerivative.DotProduct(normalUnitVec);
			return new Tuple<double[], double, double[], double[], double>(surfaceVector, m11, normalUnitVec, tangentVector, curvatureTensor);
		}

		private Tuple<double[], double, double[], double[], double> SlaveSegmentGeometry(double[,] daMatrix, double[,] da2Matrix)
		{
			var xupd = NodalXUpdated();
			var surfaceVector = Matrix.CreateFromArray(daMatrix).Multiply(xupd);
			var surfaceVectorDerivative = Matrix.CreateFromArray(da2Matrix).Multiply(xupd);
			var detm = surfaceVector.DotProduct(surfaceVector);
			var m11 = 1.0 / detm;
			var vector = new double[] { -surfaceVector[1], surfaceVector[0] };
			var scalarCoef = 1.0 / (Math.Sqrt(detm));
			var normalUnitVec = vector.Scale(scalarCoef);
			var tangentVector = surfaceVector.Scale(scalarCoef);
			var scalarCoef2 = Math.Pow(m11, 2.0);
			var curvatureTensor = scalarCoef2 * surfaceVectorDerivative.DotProduct(normalUnitVec);
			return new Tuple<double[], double, double[], double[], double>(surfaceVector, detm, normalUnitVec, tangentVector, curvatureTensor);
		}

		private double CalculatePenetration(double[,] aMatrix, double[] n)
		{
			var AT = Matrix.CreateFromArray(aMatrix).Transpose();
			var AT_n = AT.Multiply(n);
			var xupd = NodalXUpdated();
			var normalGap = xupd.DotProduct(AT_n);
			return normalGap;
		}

		private double CalculateDeltaKsi(double[] masterSlaveRelativeVector, double[] surfaceVector, double[] surfaceVectorDerivative)
		{
			var scalar1 = surfaceVector.DotProduct(masterSlaveRelativeVector);
			var scalar2 = surfaceVectorDerivative.DotProduct(masterSlaveRelativeVector) - surfaceVector.DotProduct(surfaceVector);
			var deltaKsi = -scalar1 / scalar2;
			return deltaKsi;
		}

		private double Project(double ksi1Initial, double ksi2)
		{
			if (MasterSegmentOrder == 1)
			{
				var aMatrix = CalculatePositionMatrix(ksi1Initial, ksi2).Item1;
				var m = SlaveSegmentOrder + 1;
				var slaveNMatrix = new double[2, 2 * m];
				var xUpdated = NodalXUpdated();
				var list = new List<double>();
				for (var i = 4; i < list.Count; i++)
				{
					list.Add(xUpdated[i]);
				}
				var x = list.ToArray();
				for (var i = 0; i <= 1; i++)
				{
					if (i == 0)
					{
						var countCols = 0;
						for (var j = 4; j < aMatrix.GetLength(1) - 1; j += 2)
						{
							slaveNMatrix[i, countCols] = aMatrix[i, j];
							countCols += 2;
						}
					}
					else
					{
						var countCols = 1;
						for (var j = 5; j < aMatrix.GetLength(1); j += 2)
						{
							slaveNMatrix[i, countCols] = aMatrix[i, j];
							countCols += 2;
						}
					}
				}
				var slavePositionVector = Matrix.CreateFromArray(slaveNMatrix).Multiply(x);
				var xM1 = xUpdated[0];
				var yM1 = xUpdated[1];
				var xM2 = xUpdated[2];
				var yM2 = xUpdated[3];
				var xS = slavePositionVector[0];
				var yS = slavePositionVector[1];
				var ksi = (2 * (xS * (xM2 - xM1) + yS * (yM2 - yM1)) - Math.Pow(xM2, 2) - Math.Pow(yM2, 2) + Math.Pow(xM1, 2) + Math.Pow(yM1, 2)) / (Math.Pow(xM2 - xM1, 2) + Math.Pow(yM2 - yM1, 2));
				return ksi;
			}
			else
			{
				var maxIterations = 1000;
				var tol = Math.Pow(10.0, -6.0);
				var deltaKsi = 0.0;
				var ksi = ksi1Initial;
				var xUpdated = NodalXUpdated();
				for (var i = 1; i <= maxIterations; i++)
				{
					var aMatrices = CalculatePositionMatrix(ksi, ksi2);
					var masterSlaveRelativeVector = Matrix.CreateFromArray(aMatrices.Item1).Multiply(xUpdated);
					var surfaceVector = Matrix.CreateFromArray(aMatrices.Item2).Multiply(xUpdated).Scale(-1.0);
					var surfaceVectorDerivative = Matrix.CreateFromArray(aMatrices.Item3).Multiply(xUpdated).Scale(-1.0);
					deltaKsi = CalculateDeltaKsi(masterSlaveRelativeVector, surfaceVector, surfaceVectorDerivative);
					ksi += deltaKsi;
					if (Math.Abs(deltaKsi) <= tol)
					{
						break;
					}
				}
				if (Math.Abs(deltaKsi) > tol)
				{
					throw new Exception("CPP not found in current iterations");
				}
				else
				{
					return ksi;
				}
			}
		}
		private Tuple<double[], double[]> GaussPoints()
		{
			var iP = IntegrationPoints;
			var gaussPoints = new double[iP];
			var gaussWeights = new double[iP];
			if (iP == 2)
			{
				gaussPoints[0] = -1.0 / Math.Sqrt(3);
				gaussPoints[1] = 1.0 / Math.Sqrt(3);
				gaussWeights[0] = 1.0;
				gaussWeights[1] = 1.0;
			}
			else if (iP == 3)
			{
				gaussPoints[0] = -0.77459;
				gaussPoints[1] = 0.0;
				gaussPoints[2] = 0.77459;
				gaussWeights[0] = 0.55555;
				gaussWeights[1] = 0.88888;
				gaussWeights[2] = 0.55555;
			}
			else if (iP == 4)
			{
				gaussPoints[0] = -0.86113;
				gaussPoints[1] = -0.33998;
				gaussPoints[2] = 0.33998;
				gaussPoints[3] = 0.86113;
				gaussWeights[0] = 0.34785;
				gaussWeights[1] = 0.65214;
				gaussWeights[2] = 0.65214;
				gaussWeights[3] = 0.34785;
			}
			else if (iP == 5)
			{
				gaussPoints[0] = -0.90617;
				gaussPoints[1] = -0.53846;
				gaussPoints[2] = 0.0;
				gaussPoints[3] = 0.53846;
				gaussPoints[4] = 0.90617;
				gaussWeights[0] = 0.23692;
				gaussWeights[1] = 0.47862;
				gaussWeights[2] = 0.56888;
				gaussWeights[3] = 0.47862;
				gaussWeights[4] = 0.23692;
			}
			else if (iP == 6)
			{
				gaussPoints[0] = -0.93246;
				gaussPoints[1] = -0.66120;
				gaussPoints[2] = -0.23861;
				gaussPoints[3] = 0.23861;
				gaussPoints[4] = 0.66120;
				gaussPoints[5] = 0.93246;
				gaussWeights[0] = 0.17132;
				gaussWeights[1] = 0.36076;
				gaussWeights[2] = 0.46791;
				gaussWeights[3] = 0.46791;
				gaussWeights[4] = 0.36076;
				gaussWeights[5] = 0.17132;
			}
			return new Tuple<double[], double[]>(gaussPoints, gaussWeights);
		}

		#region IElementType Members

		public int ID { get; set; }
		public int SubdomainID { get; set; }
		public IReadOnlyList<INode> Nodes { get; }

		public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofs;

		public IReadOnlyList<INode> GetNodesForMatrixAssembly() => Nodes;

		private Matrix CalculateMainStiffnessPart(double ksi1, double ksi2, double[] n)
		{
			var positionMatrices = CalculatePositionMatrix(ksi1, ksi2);
			var A = positionMatrices.Item1;
			var nxn = n.TensorProduct(n);
			var nxn_A = nxn * Matrix.CreateFromArray(A);
			var AT_nxn_A = Matrix.CreateFromArray(A).Transpose() * nxn_A;
			var mainStiffnessMatrix = AT_nxn_A.Scale(penaltyFactor);
			return mainStiffnessMatrix;
		}

		private Matrix CalculateRotationalStiffnessPart(double[,] A, double[,] dA, double[] n, double ksi3, double m11, double[] dRho)
		{
			var coef = penaltyFactor * ksi3 * m11;
			var n_x_dRho = n.TensorProduct(dRho);
			var dRho_x_n = dRho.TensorProduct(n);
			var firstTerm = Matrix.CreateFromArray(dA).Transpose() * n_x_dRho * Matrix.CreateFromArray(A);
			var secondTerm = Matrix.CreateFromArray(A).Transpose() * dRho_x_n * Matrix.CreateFromArray(dA);
			var rotationalPart = (firstTerm + secondTerm).Scale(coef);
			return rotationalPart;
		}

		private Matrix CalculateCurvatureStiffnessPart(double[,] A, double ksi3, double m11, double[] dRho, double h11)
		{
			var coef = penaltyFactor * ksi3 * m11 * h11;
			var dRho_x_dRho = dRho.TensorProduct(dRho);
			var mat = Matrix.CreateFromArray(A).Transpose() * dRho_x_dRho * Matrix.CreateFromArray(A);
			var curvaturePart = mat.Scale(coef);
			return curvaturePart;
		}

		public IMatrix StiffnessMatrix()
		{
			var globalStifnessMatrix = Matrix.CreateFromArray(new double[2 * Nodes.Count, 2 * Nodes.Count]);
			var gPArray = GaussPoints().Item1;
			var gWArray = GaussPoints().Item2;
			for (var i = 0; i < IntegrationPoints; i++)
			{
				var ksi2 = gPArray[i];
				var gW = gWArray[i];
				var ksi1 = Project(0.0, ksi2);
				if (Math.Abs(ksi1) <= 1.05)
				{
					var positionMatrices = CalculatePositionMatrix(ksi1, ksi2);
					var aMatrix = positionMatrices.Item1;
					var daMatrix = positionMatrices.Item2;
					var da2Matrix = positionMatrices.Item3;
					var masterSurfaceCharacteristics = MasterSegmentGeometry(daMatrix, da2Matrix);
					var m11 = masterSurfaceCharacteristics.Item2;
					var dRho = masterSurfaceCharacteristics.Item1;
					var n = masterSurfaceCharacteristics.Item3;
					var h11 = masterSurfaceCharacteristics.Item5;
					var ksi3 = CalculatePenetration(aMatrix, n);
					if (ksi3 <= 0)
					{
						var slaveMetricTensor = SlaveSegmentGeometry(positionMatrices.Item4, positionMatrices.Item5).Item2;
						var mainPart = CalculateMainStiffnessPart(ksi1, ksi2, n);
						var rotationalPart = CalculateRotationalStiffnessPart(aMatrix, daMatrix, n, ksi3, m11, dRho);
						var curvaturePart = CalculateCurvatureStiffnessPart(aMatrix, ksi3, m11, dRho, h11);
						var scalar = Math.Pow(slaveMetricTensor, 0.5) * gW;
						var StifnessMatrix = (mainPart + rotationalPart + curvaturePart).Scale(scalar);
						globalStifnessMatrix += StifnessMatrix;
					}
				}
			}
			return dofEnumerator.GetTransformedMatrix(globalStifnessMatrix);
		}

		public IMatrix PhysicsMatrix() => StiffnessMatrix();

		public IMatrix MassMatrix()
		{
			double[,] massMatrix = new double[2 * Nodes.Count, 2 * Nodes.Count];
			return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(massMatrix));
		}

		public IMatrix DampingMatrix()
		{
			double[,] dampingMatrix = new double[2 * Nodes.Count, 2 * Nodes.Count];
			return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(dampingMatrix));
		}

		private double[] CreateInternalGlobalForcesVector()
		{
			var internalGlobalForcesVector = new double[2 * Nodes.Count];
			var gPArray = GaussPoints().Item1;
			var gWArray = GaussPoints().Item2;
			for (var i = 0; i < IntegrationPoints; i++)
			{
				var ksi2 = gPArray[i];
				var gW = gWArray[i];
				var ksi1 = Project(0.0, ksi2);
				if (Math.Abs(ksi1) <= 1.05)
				{
					var positionMatrices = CalculatePositionMatrix(ksi1, ksi2);
					var aMatrix = positionMatrices.Item1;
					var daMatrix = positionMatrices.Item2;
					var da2Matrix = positionMatrices.Item3;
					var surfaceCharacteristics = MasterSegmentGeometry(daMatrix, da2Matrix);
					var n = surfaceCharacteristics.Item3;
					var ksi3 = CalculatePenetration(aMatrix, n);
					if (ksi3 <= 0)
					{
						var slaveMetricTensor = SlaveSegmentGeometry(positionMatrices.Item4, positionMatrices.Item5).Item2;
						var scalar = Math.Pow(slaveMetricTensor, 0.5) * gW;
						var AT = Matrix.CreateFromArray(aMatrix).Transpose();
						var AT_n = AT.Multiply(n);
						var internalForcesVectorGaussPoint = AT_n.Scale(penaltyFactor * ksi3 * scalar);
						internalGlobalForcesVector = internalGlobalForcesVector.Add(internalForcesVectorGaussPoint);
					}
				}
			}
			return internalGlobalForcesVector;
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
