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
using MGroup.LinearAlgebra;

namespace MGroup.FEM.Structural.Line
{
	public class ContactSurfaceToSurface3DSerendipitySecondOrder : IStructuralElementType
	{
		private readonly IDofType[][] dofs;
		private readonly double penaltyFactor;
		private readonly int MasterSurfaceOrder;
		private readonly int SlaveSurfaceOrder;
		private readonly int IntegrationPointsPerNaturalAxis;

		private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
		private double[] DisplacementVector { get; set; }
		private double ContactArea { get; }
		//private List<int> activeElements = new List<int>();
		public ContactSurfaceToSurface3DSerendipitySecondOrder(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplier, double contactArea)
		{
			if (nodes.Count != 8)
			{
				throw new ArgumentException("This Constructor can be used only for linear Surfaces");
			}
			else
			{
				this.penaltyFactor = penaltyFactorMultiplier * youngModulus;
				this.DisplacementVector = new double[3 * nodes.Count];
				this.ContactArea = contactArea;
				this.Nodes = nodes;
				dofs = new IDofType[nodes.Count][];
				for (var i = 0; i < nodes.Count; i++)
				{
					dofs[i] = new IDofType[]
					{
						StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
					};
				}
				this.MasterSurfaceOrder = 1;
				this.SlaveSurfaceOrder = 1;
				this.IntegrationPointsPerNaturalAxis = 2;
			}
		}
		public ContactSurfaceToSurface3DSerendipitySecondOrder(IReadOnlyList<INode> nodes, double penaltyFactor, double contactArea)
		{
			if (nodes.Count != 8)
			{
				throw new ArgumentException("This Constructor can be used only for linear Surfaces");
			}
			this.penaltyFactor = penaltyFactor;
			this.DisplacementVector = new double[3 * nodes.Count];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
			dofs = new IDofType[nodes.Count][];
			for (var i = 0; i < nodes.Count; i++)
			{
				dofs[i] = new IDofType[]
				{
						StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
				};
			}
			this.MasterSurfaceOrder = 1;
			this.SlaveSurfaceOrder = 1;
			this.IntegrationPointsPerNaturalAxis = 2;
		}
		public ContactSurfaceToSurface3DSerendipitySecondOrder(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplier, double contactArea,
			int masterSurfaceOrder, int slaveSurfaceOrder)
		{
			if ((masterSurfaceOrder != 1 && masterSurfaceOrder != 2) || (slaveSurfaceOrder != 1 && slaveSurfaceOrder != 2))
			{
				throw new ArgumentException("Only linear & quadratic Surfaces can be defined");
			}
			if (((masterSurfaceOrder == 1 && slaveSurfaceOrder == 1) && nodes.Count != 8) ||
				((masterSurfaceOrder == 1 && slaveSurfaceOrder == 2) && nodes.Count != 12) ||
				((masterSurfaceOrder == 2 && slaveSurfaceOrder == 1) && nodes.Count != 12) ||
				((masterSurfaceOrder == 2 && slaveSurfaceOrder == 2) && nodes.Count != 16))
			{
				throw new ArgumentException("Inconsistent input regarding the nodes & the order of the Surfaces");
			}
			this.penaltyFactor = penaltyFactorMultiplier * youngModulus;
			this.DisplacementVector = new double[3 * nodes.Count];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
			dofs = new IDofType[nodes.Count][];
			for (var i = 0; i < nodes.Count; i++)
			{
				dofs[i] = new IDofType[]
				{
							StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
				};
			}
			this.MasterSurfaceOrder = masterSurfaceOrder;
			this.SlaveSurfaceOrder = slaveSurfaceOrder;
			if (slaveSurfaceOrder == 1)
			{
				this.IntegrationPointsPerNaturalAxis = 2;
			}
			else
			{
				this.IntegrationPointsPerNaturalAxis = 3;
			}
		}
		public ContactSurfaceToSurface3DSerendipitySecondOrder(IReadOnlyList<INode> nodes, double penaltyFactor, double contactArea, int masterSurfaceOrder, int slaveSurfaceOrder)
		{
			if ((masterSurfaceOrder != 1 && masterSurfaceOrder != 2) || (slaveSurfaceOrder != 1 && slaveSurfaceOrder != 2))
			{
				throw new ArgumentException("Only linear & quadratic Surfaces can be defined");
			}
			if (((masterSurfaceOrder == 1 && slaveSurfaceOrder == 1) && nodes.Count != 8) ||
				((masterSurfaceOrder == 1 && slaveSurfaceOrder == 2) && nodes.Count != 12) ||
				((masterSurfaceOrder == 2 && slaveSurfaceOrder == 1) && nodes.Count != 12) ||
				((masterSurfaceOrder == 2 && slaveSurfaceOrder == 2) && nodes.Count != 16))
			{
				throw new ArgumentException("Inconsistent input regarding the nodes & the order of the Surfaces");
			}
			this.penaltyFactor = penaltyFactor;
			this.DisplacementVector = new double[3 * nodes.Count];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
			dofs = new IDofType[nodes.Count][];
			for (var i = 0; i < nodes.Count; i++)
			{
				dofs[i] = new IDofType[]
				{
							StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
				};
			}
			this.MasterSurfaceOrder = masterSurfaceOrder;
			this.SlaveSurfaceOrder = slaveSurfaceOrder;
			if (slaveSurfaceOrder == 1)
			{
				this.IntegrationPointsPerNaturalAxis = 2;
			}
			else
			{
				this.IntegrationPointsPerNaturalAxis = 3;
			}
		}
		public ContactSurfaceToSurface3DSerendipitySecondOrder(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplier, double contactArea,
			int masterSurfaceOrder, int slaveSurfaceOrder, int integrationPointsPerNaturalAxis)
		{
			if ((masterSurfaceOrder != 1 && masterSurfaceOrder != 2) || (slaveSurfaceOrder != 1 && slaveSurfaceOrder != 2))
			{
				throw new ArgumentException("Only linear & quadratic Surfaces can be defined");
			}
			if (((masterSurfaceOrder == 1 && slaveSurfaceOrder == 1) && nodes.Count != 8) ||
				((masterSurfaceOrder == 1 && slaveSurfaceOrder == 2) && nodes.Count != 12) ||
				((masterSurfaceOrder == 2 && slaveSurfaceOrder == 1) && nodes.Count != 12) ||
				((masterSurfaceOrder == 2 && slaveSurfaceOrder == 2) && nodes.Count != 16))
			{
				throw new ArgumentException("Inconsistent input regarding the nodes & the order of the Surfaces");
			}
			if (integrationPointsPerNaturalAxis < 2 || integrationPointsPerNaturalAxis > 10)
			{
				throw new ArgumentException("Between [2,10] Gauss points per Natural Axis can be defined");
			}
			this.penaltyFactor = penaltyFactorMultiplier * youngModulus;
			this.DisplacementVector = new double[3 * nodes.Count];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
			dofs = new IDofType[nodes.Count][];
			for (var i = 0; i < nodes.Count; i++)
			{
				dofs[i] = new IDofType[]
				{
							StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
				};
			}
			this.MasterSurfaceOrder = masterSurfaceOrder;
			this.SlaveSurfaceOrder = slaveSurfaceOrder;
			this.IntegrationPointsPerNaturalAxis = integrationPointsPerNaturalAxis;
		}
		public ContactSurfaceToSurface3DSerendipitySecondOrder(IReadOnlyList<INode> nodes, double penaltyFactor, double contactArea, int masterSurfaceOrder, int slaveSurfaceOrder, int integrationPointsPerNaturalAxis)
		{
			if ((masterSurfaceOrder != 1 && masterSurfaceOrder != 2) || (slaveSurfaceOrder != 1 && slaveSurfaceOrder != 2))
			{
				throw new ArgumentException("Only linear & quadratic Surfaces can be defined");
			}
			if (((masterSurfaceOrder == 1 && slaveSurfaceOrder == 1) && nodes.Count != 8) ||
				((masterSurfaceOrder == 1 && slaveSurfaceOrder == 2) && nodes.Count != 12) ||
				((masterSurfaceOrder == 2 && slaveSurfaceOrder == 1) && nodes.Count != 12) ||
				((masterSurfaceOrder == 2 && slaveSurfaceOrder == 2) && nodes.Count != 16))
			{
				throw new ArgumentException("Inconsistent input regarding the nodes & the order of the Surfaces");
			}
			if (integrationPointsPerNaturalAxis < 2 || integrationPointsPerNaturalAxis > 10)
			{
				throw new ArgumentException("Between [2,10] Gauss points per Natural Axis can be defined");
			}
			this.penaltyFactor = penaltyFactor;
			this.DisplacementVector = new double[3 * nodes.Count];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
			dofs = new IDofType[nodes.Count][];
			for (var i = 0; i < nodes.Count; i++)
			{
				dofs[i] = new IDofType[]
				{
							StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ
				};
			}
			this.MasterSurfaceOrder = masterSurfaceOrder;
			this.SlaveSurfaceOrder = slaveSurfaceOrder;
			this.IntegrationPointsPerNaturalAxis = integrationPointsPerNaturalAxis;
		}
		public ContactSurfaceToSurface3DSerendipitySecondOrder(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplier,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, youngModulus, penaltyFactorMultiplier, 1d)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public ContactSurfaceToSurface3DSerendipitySecondOrder(IReadOnlyList<INode> nodes, double penaltyFactor, IElementDofEnumerator dofEnumerator)
			: this(nodes, penaltyFactor, 1.0)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public ContactSurfaceToSurface3DSerendipitySecondOrder(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplier, int masterSurfaceOrder, int slaveSurfaceOrder,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, youngModulus, penaltyFactorMultiplier, 1d, masterSurfaceOrder, slaveSurfaceOrder)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public ContactSurfaceToSurface3DSerendipitySecondOrder(IReadOnlyList<INode> nodes, double penaltyFactor, int masterSurfaceOrder, int slaveSurfaceOrder,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, penaltyFactor, 1.0, masterSurfaceOrder, slaveSurfaceOrder)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public ContactSurfaceToSurface3DSerendipitySecondOrder(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplier, int masterSurfaceOrder, int slaveSurfaceOrder, int integrationPointsPerNaturalAxis,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, youngModulus, penaltyFactorMultiplier, 1d, masterSurfaceOrder, slaveSurfaceOrder, integrationPointsPerNaturalAxis)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public ContactSurfaceToSurface3DSerendipitySecondOrder(IReadOnlyList<INode> nodes, double penaltyFactor, int masterSurfaceOrder, int slaveSurfaceOrder, int integrationPointsPerNaturalAxis,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, penaltyFactor, 1.0, masterSurfaceOrder, slaveSurfaceOrder, integrationPointsPerNaturalAxis)
		{
			this.dofEnumerator = dofEnumerator;
		}

		public CellType CellType { get; } = CellType.Unknown;

		public IElementDofEnumerator DofEnumerator
		{
			get { return dofEnumerator; }
			set { dofEnumerator = value; }
		}
		private double[] NodalXUpdated()
		{
			var x = new double[3 * Nodes.Count];
			for (var i = 0; i < Nodes.Count; i++)
			{
				x[3 * i] = Nodes[i].X + DisplacementVector[3 * i];
				x[3 * i + 1] = Nodes[i].Y + DisplacementVector[3 * i + 1];
				x[3 * i + 2] = Nodes[i].Z + DisplacementVector[3 * i + 2];
			}
			return x;
		}
		private Tuple<double[,], Tuple<double[,], double[,], double[,], double[,], double[,]>, Tuple<double[,], double[,], double[,], double[,]>> CalculatePositionMatrix(double ksi1, double ksi2, double ksi3, double ksi4)
		{
			if (SlaveSurfaceOrder == 1 && MasterSurfaceOrder == 1)
			{
				var N1 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 - ksi2);
				var N2 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 - ksi2);
				var N3 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 + ksi2);
				var N4 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 + ksi2);
				var N5 = 1.0 / 4.0 * (1.0 - ksi3) * (1.0 - ksi4);
				var N6 = 1.0 / 4.0 * (1.0 + ksi3) * (1.0 - ksi4);
				var N7 = 1.0 / 4.0 * (1.0 + ksi3) * (1.0 + ksi4);
				var N8 = 1.0 / 4.0 * (1.0 - ksi3) * (1.0 + ksi4);

				var dN11 = -1.0 / 4.0 * (1.0 - ksi2);
				var dN21 = 1.0 / 4.0 * (1.0 - ksi2);
				var dN31 = 1.0 / 4.0 * (1.0 + ksi2);
				var dN41 = -1.0 / 4.0 * (1.0 + ksi2);

				var dN12 = -1.0 / 4.0 * (1.0 - ksi1);
				var dN22 = -1.0 / 4.0 * (1.0 + ksi1);
				var dN32 = 1.0 / 4.0 * (1.0 + ksi1);
				var dN42 = 1.0 / 4.0 * (1.0 - ksi1);

				var dN53 = -1.0 / 4.0 * (1.0 - ksi4);
				var dN63 = 1.0 / 4.0 * (1.0 - ksi4);
				var dN73 = 1.0 / 4.0 * (1.0 + ksi4);
				var dN83 = -1.0 / 4.0 * (1.0 + ksi4);

				var dN54 = -1.0 / 4.0 * (1.0 - ksi3);
				var dN64 = -1.0 / 4.0 * (1.0 + ksi3);
				var dN74 = 1.0 / 4.0 * (1.0 + ksi3);
				var dN84 = 1.0 / 4.0 * (1.0 - ksi3);

				var dN112 = 1.0 / 4.0;
				var dN212 = -1.0 / 4.0;
				var dN312 = 1.0 / 4.0;
				var dN412 = -1.0 / 4.0;

				var aMatrix = new double[,]
					{
					{ -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, N5 ,0.0, 0.0 , N6, 0.0 ,0.0 , N7, 0.0, 0.0, N8, 0.0, 0.0 },
					{ 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, N5 ,0.0, 0.0 , N6, 0.0 ,0.0 , N7, 0.0, 0.0, N8, 0.0 },
					{ 0.0, 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, N5 ,0.0, 0.0 , N6, 0.0 ,0.0 , N7, 0.0, 0.0, N8 }
					};

				var da1Matrix = new double[,]
					{
					{ -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0 }
					};

				var da11Matrix = new double[,]
					{
					{ 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0 }
					};

				var da12Matrix = new double[,]
					{
					{ -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0 }
					};

				var da2Matrix = new double[,]
					{
					{ -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0 }
					};

				var da22Matrix = new double[,]
					{
					{ 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0 }
					};

				var da3Matrix = new double[,]
					{
					{ 0.0, 0.0, 0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN53 ,0.0, 0.0 , dN63, 0.0 ,0.0 , dN73, 0.0, 0.0, dN83, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN53 ,0.0, 0.0 , dN63, 0.0 ,0.0 , dN73, 0.0, 0.0, dN83, 0.0 },
					{ 0.0, 0.0, 0.0 ,0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN53 ,0.0, 0.0 , dN63, 0.0 ,0.0 , dN73, 0.0, 0.0, dN83 }
					};

				var da33Matrix = new double[,]
					{
					{ 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0 }
					};

				var da4Matrix = new double[,]
					{
					{ 0.0, 0.0, 0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN54 ,0.0, 0.0 , dN64, 0.0 ,0.0 , dN74, 0.0, 0.0, dN84, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN54 ,0.0, 0.0 , dN64, 0.0 ,0.0 , dN74, 0.0, 0.0, dN84, 0.0 },
					{ 0.0, 0.0, 0.0 ,0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN54 ,0.0, 0.0 , dN64, 0.0 ,0.0 , dN74, 0.0, 0.0, dN84 }
					};

				var da44Matrix = new double[,]
					{
					{ 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0 ,0.0, 0.0 ,0.0 ,0.0 ,0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0 }
					};

				var T1 = new Tuple<double[,], double[,], double[,], double[,], double[,]>(da1Matrix, da11Matrix, da2Matrix, da22Matrix, da12Matrix);
				var T2 = new Tuple<double[,], double[,], double[,], double[,]>(da3Matrix, da33Matrix, da4Matrix, da44Matrix);
				return new Tuple<double[,], Tuple<double[,], double[,], double[,], double[,], double[,]>, Tuple<double[,], double[,], double[,], double[,]>>(aMatrix, T1, T2);
			}
			else if (MasterSurfaceOrder == 1 && SlaveSurfaceOrder == 2)
			{
				var N1 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 - ksi2);
				var N2 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 - ksi2);
				var N3 = 1.0 / 4.0 * (1.0 + ksi1) * (1.0 + ksi2);
				var N4 = 1.0 / 4.0 * (1.0 - ksi1) * (1.0 + ksi2);

				var N5 = 1.0 / 4.0 * (1d + ksi3) * (1d + ksi4) - 1.0 / 4.0 * (1d + ksi3) * (1 - ksi4 * ksi4) - 1.0 / 4.0 * (1d - ksi3 * ksi3) * (1d + ksi4 );
				var N6 = 1.0 / 4.0 * (1d - ksi3) * (1d + ksi4) - 1.0 / 4.0 * (1d - ksi3 * ksi3) * (1 + ksi4) - 1.0 / 4.0 * (1d - ksi3) * (1d - ksi4 * ksi4);
				var N7 = 1.0 / 4.0 * (1d - ksi3) * (1d - ksi4) - 1.0 / 4.0 * (1d - ksi3 * ksi3) * (1 - ksi4) - 1.0 / 4.0 * (1d - ksi3) * (1d - ksi4 * ksi4);
				var N8 = 1.0 / 4.0 * (1d + ksi3) * (1d - ksi4) - 1.0 / 4.0 * (1d - ksi3 * ksi3) * (1 - ksi4) - 1.0 / 4.0 * (1d + ksi3) * (1d - ksi4 * ksi4);
				var N9 = 1.0 / 2.0 * (1d - ksi3 * ksi3) * (1d + ksi4);
				var N10 = 1.0 / 2.0 * (1d - ksi3) * (1d - ksi4 * ksi4);
				var N11 = 1.0 / 2.0 * (1d - ksi3 * ksi3) * (1d - ksi4);
				var N12 = 1.0 / 2.0 * (1d + ksi3) * (1d - ksi4 * ksi4);

				var dN11 = -1.0 / 4.0 * (1.0 - ksi2);
				var dN21 = 1.0 / 4.0 * (1.0 - ksi2);
				var dN31 = 1.0 / 4.0 * (1.0 + ksi2);
				var dN41 = -1.0 / 4.0 * (1.0 + ksi2);

				var dN12 = -1.0 / 4.0 * (1.0 - ksi1);
				var dN22 = -1.0 / 4.0 * (1.0 + ksi1);
				var dN32 = 1.0 / 4.0 * (1.0 + ksi1);
				var dN42 = 1.0 / 4.0 * (1.0 - ksi1);

				var dN112 = 1.0 / 4.0;
				var dN212 = -1.0 / 4.0;
				var dN312 = 1.0 / 4.0;
				var dN412 = -1.0 / 4.0;

				var dN53 = 1.0 / 4.0 * (2d * ksi3 + ksi4) * (1d + ksi4);
				var dN63 = 1.0 / 4.0 * (2d * ksi3 - ksi4) * (1d + ksi4);
				var dN73 = 1.0 / 4.0 * (2d * ksi3 + ksi4) * (1d - ksi4);
				var dN83 = 1.0 / 4.0 * (2d * ksi3 - ksi4) * (1d - ksi4);
				var dN93 = -ksi3 * (1d + ksi4);
				var dN103 = -1.0 / 2.0 * (1d - ksi4* ksi4);
				var dN113 = -ksi3 * (1d - ksi4);
				var dN123 = 1.0 / 2.0 * (1d - ksi4 * ksi4);

				var dN533 = 1.0 / 2.0 * (1d + ksi4);
				var dN633 = 1.0 / 2.0 * (1d + ksi4);
				var dN733 = 1.0 / 2.0 * (1d - ksi4);
				var dN833 = 1.0 / 2.0 * (1d - ksi4);
				var dN933 = -(1d + ksi4);
				var dN1033 = 0d;
				var dN1133 = (ksi4 - 1d);
				var dN1233 = 0d;

				var dN54 = 1.0 / 4.0 * (ksi3 + 2d * ksi4) * (1d + ksi3);
				var dN64 = 1.0 / 4.0 * (-ksi3 + 2d * ksi4) * (1d - ksi3);
				var dN74 = 1.0 / 4.0 * (ksi3 + 2d * ksi4) * (1d - ksi3);
				var dN84 = 1.0 / 4.0 * (-ksi3 + 2d * ksi4) * (1d + ksi3);
				var dN94 = 1.0 / 2.0 * (1d - Math.Pow(ksi3, 2));
				var dN104 = -(1d - ksi3) * ksi4;
				var dN114 = -1.0 / 2.0 * (1d - Math.Pow(ksi3, 2));
				var dN124 = -(1d + ksi3) * ksi4;

				var dN544 = 1.0 / 2.0 * (1d + ksi3);
				var dN644 = 1.0 / 2.0 * (1d - ksi3);
				var dN744 = 1.0 / 2.0 * (1d - ksi3);
				var dN844 = 1.0 / 2.0 * (1d + ksi3);
				var dN944 = 0d;
				var dN1044 = ksi3 - 1d;
				var dN1144 = 0d;
				var dN1244 = -1d - ksi3;

				var aMatrix = new double[,]
					{
					{ -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, N5 ,0.0, 0.0 , N6, 0.0 ,0.0 , N7, 0.0, 0.0, N8, 0.0, 0.0, N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12, 0.0, 0.0},
					{ 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, N5 ,0.0, 0.0 , N6, 0.0 ,0.0 , N7, 0.0, 0.0, N8, 0.0,0.0, N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12, 0.0 },
					{ 0.0, 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, N5 ,0.0, 0.0 , N6, 0.0 ,0.0 , N7, 0.0, 0.0, N8, 0.0,0.0, N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12 }
					};

				var da1Matrix = new double[,]
					{
					{ -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
					};

				var da11Matrix = new double[,]
					{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
					};

				var da12Matrix = new double[,]
					{
					{ -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
					 };
				var da2Matrix = new double[,]
					{
					{ -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
					};

				var da22Matrix = new double[,]
					{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
					};

				var da3Matrix = new double[,]
					{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN53 ,0.0, 0.0 , dN63, 0.0 ,0.0 , dN73, 0.0, 0.0, dN83, 0.0, 0.0, dN93, 0.0, 0.0, dN103, 0.0, 0.0, dN113, 0.0, 0.0, dN123, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN53 ,0.0, 0.0 , dN63, 0.0 ,0.0 , dN73, 0.0, 0.0, dN83, 0.0,0.0, dN93, 0.0, 0.0, dN103, 0.0, 0.0, dN113, 0.0, 0.0, dN123, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN53 ,0.0, 0.0 , dN63, 0.0 ,0.0 , dN73, 0.0, 0.0, dN83, 0.0,0.0, dN93, 0.0, 0.0, dN103, 0.0, 0.0, dN113, 0.0, 0.0, dN123 }
					};

				var da33Matrix = new double[,]
					{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN533 ,0.0, 0.0 , dN633, 0.0 ,0.0 , dN733, 0.0, 0.0, dN833, 0.0, 0.0, dN933, 0.0, 0.0, dN1033, 0.0, 0.0, dN1133, 0.0, 0.0, dN1233, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN533 ,0.0, 0.0 , dN633, 0.0 ,0.0 , dN733, 0.0, 0.0, dN833, 0.0,0.0, dN933, 0.0, 0.0, dN1033, 0.0, 0.0, dN1133, 0.0, 0.0, dN1233, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN533 ,0.0, 0.0 , dN633, 0.0 ,0.0 , dN733, 0.0, 0.0, dN833, 0.0,0.0, dN933, 0.0, 0.0, dN1033, 0.0, 0.0, dN1133, 0.0, 0.0, dN1233 }
					};

				var da4Matrix = new double[,]
					{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN54 ,0.0, 0.0 , dN64, 0.0 ,0.0 , dN74, 0.0, 0.0, dN84, 0.0, 0.0, dN94, 0.0, 0.0, dN104, 0.0, 0.0, dN114, 0.0, 0.0, dN124, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN54 ,0.0, 0.0 , dN64, 0.0 ,0.0 , dN74, 0.0, 0.0, dN84, 0.0,0.0, dN94, 0.0, 0.0, dN104, 0.0, 0.0, dN114, 0.0, 0.0, dN124, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN54 ,0.0, 0.0 , dN64, 0.0 ,0.0 , dN74, 0.0, 0.0, dN84, 0.0,0.0, dN94, 0.0, 0.0, dN104, 0.0, 0.0, dN114, 0.0, 0.0, dN124 }
					};

				var da44Matrix = new double[,]
					{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN544 ,0.0, 0.0 , dN644, 0.0 ,0.0 , dN744, 0.0, 0.0, dN844, 0.0, 0.0, dN944, 0.0, 0.0, dN1044, 0.0, 0.0, dN1144, 0.0, 0.0, dN1244, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN544 ,0.0, 0.0 , dN644, 0.0 ,0.0 , dN744, 0.0, 0.0, dN844, 0.0,0.0, dN944, 0.0, 0.0, dN1044, 0.0, 0.0, dN1144, 0.0, 0.0, dN1244, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN544 ,0.0, 0.0 , dN644, 0.0 ,0.0 , dN744, 0.0, 0.0, dN844, 0.0,0.0, dN944, 0.0, 0.0, dN1044, 0.0, 0.0, dN1144, 0.0, 0.0, dN1244 }
					};

				var T1 = new Tuple<double[,], double[,], double[,], double[,], double[,]>(da1Matrix, da11Matrix, da2Matrix, da22Matrix, da12Matrix);
				var T2 = new Tuple<double[,], double[,], double[,], double[,]>(da3Matrix, da33Matrix, da4Matrix, da44Matrix);
				return new Tuple<double[,], Tuple<double[,], double[,], double[,], double[,], double[,]>, Tuple<double[,], double[,], double[,], double[,]>>(aMatrix, T1, T2);
			}
			else if (MasterSurfaceOrder == 2 && SlaveSurfaceOrder == 1)
			{
				
				var N1 = 1.0 / 4.0 * (1d + ksi1) * (1d + ksi2) - 1.0 / 4.0 * (1d + ksi1) * (1 - ksi2 * ksi2) - 1.0 / 4.0 * (1d - ksi1 * ksi1) * (1d + ksi2);
				var N2 = 1.0 / 4.0 * (1d - ksi1) * (1d + ksi2) - 1.0 / 4.0 * (1d - ksi1 * ksi1) * (1 + ksi2) - 1.0 / 4.0 * (1d - ksi1) * (1d - ksi2 * ksi2);
				var N3 = 1.0 / 4.0 * (1d - ksi1) * (1d - ksi2) - 1.0 / 4.0 * (1d - ksi1 * ksi1) * (1 - ksi2) - 1.0 / 4.0 * (1d - ksi1) * (1d - ksi2 * ksi2);
				var N4 = 1.0 / 4.0 * (1d + ksi1) * (1d - ksi2) - 1.0 / 4.0 * (1d - ksi1 * ksi1) * (1 - ksi2) - 1.0 / 4.0 * (1d + ksi1) * (1d - ksi2 * ksi2);
				var N5 = 1.0 / 2.0 * (1d - ksi1 * ksi1) * (1d + ksi2);
				var N6 = 1.0 / 2.0 * (1d - ksi1) * (1d - ksi2 * ksi2);
				var N7 = 1.0 / 2.0 * (1d - ksi1 * ksi1) * (1d - ksi2);
				var N8 = 1.0 / 2.0 * (1d + ksi1) * (1d - ksi2 * ksi2);

				var dN11 = 1.0 / 4.0 * (2d * ksi1 + ksi2) * (1d + ksi2);
				var dN21 = 1.0 / 4.0 * (2d * ksi1 - ksi2) * (1d + ksi2);
				var dN31 = 1.0 / 4.0 * (2d * ksi1 + ksi2) * (1d - ksi2);
				var dN41 = 1.0 / 4.0 * (2d * ksi1 - ksi2) * (1d - ksi2);
				var dN51 = -ksi1 * (1d + ksi2);
				var dN61 = -1.0 / 2.0 * (1d - ksi2 * ksi2);
				var dN71 = -ksi1 * (1d - ksi2);
				var dN81 = 1.0 / 2.0 * (1d - ksi2 * ksi2);

				var dN111 = 1.0 / 2.0 * (1d + ksi2);
				var dN211 = 1.0 / 2.0 * (1d + ksi2);
				var dN311 = 1.0 / 2.0 * (1d - ksi2);
				var dN411 = 1.0 / 2.0 * (1d - ksi2);
				var dN511 = -(1d + ksi2);
				var dN611 = 0d;
				var dN711 = (ksi2 - 1d);
				var dN811 = 0d;

				var dN12 = 1.0 / 4.0 * (ksi1 + 2d * ksi2) * (1d + ksi1);
				var dN22 = 1.0 / 4.0 * (-ksi1 + 2d * ksi2) * (1d - ksi1);
				var dN32 = 1.0 / 4.0 * (ksi1 + 2d * ksi2) * (1d - ksi1);
				var dN42 = 1.0 / 4.0 * (-ksi1 + 2d * ksi2) * (1d + ksi1);
				var dN52 = 1.0 / 2.0 * (1d - Math.Pow(ksi1, 2));
				var dN62 = -(1d - ksi1) * ksi2;
				var dN72 = -1.0 / 2.0 * (1d - Math.Pow(ksi1, 2));
				var dN82 = -(1d + ksi1) * ksi2;

				var dN122 = 1.0 / 2.0 * (1d + ksi1);
				var dN222 = 1.0 / 2.0 * (1d - ksi1);
				var dN322 = 1.0 / 2.0 * (1d - ksi1);
				var dN422 = 1.0 / 2.0 * (1d + ksi1);
				var dN522 = 0d;
				var dN622 = ksi1 - 1d;
				var dN722 = 0d;
				var dN822 = -1d - ksi1;

				var dN112 = 1.0 / 4.0 * (2d * ksi1 + 2d * ksi2 + 1d);
				var dN212 = 1.0 / 4.0 * (2d * ksi1 - 2d * ksi2 - 1d);
				var dN312 = 1.0 / 4.0 * (-2d * ksi1 - 2d * ksi2 + 1d);
				var dN412 = 1.0 / 4.0 * (-2d * ksi1 + 2d * ksi2 - 1d);
				var dN512 = -ksi1;
				var dN612 = ksi2;
				var dN712 = ksi1;
				var dN812 = -ksi2;


				var N9 = 1.0 / 4.0 * (1.0 - ksi3) * (1.0 - ksi4);
				var N10 = 1.0 / 4.0 * (1.0 + ksi3) * (1.0 - ksi4);
				var N11 = 1.0 / 4.0 * (1.0 + ksi3) * (1.0 + ksi4);
				var N12 = 1.0 / 4.0 * (1.0 - ksi3) * (1.0 + ksi4);

				var dN93 = -1.0 / 4.0 * (1.0 - ksi4);
				var dN103 = 1.0 / 4.0 * (1.0 - ksi4);
				var dN113 = 1.0 / 4.0 * (1.0 + ksi4);
				var dN123 = -1.0 / 4.0 * (1.0 + ksi4);

				var dN94 = -1.0 / 4.0 * (1.0 - ksi3);
				var dN104 = -1.0 / 4.0 * (1.0 + ksi3);
				var dN114 = 1.0 / 4.0 * (1.0 + ksi3);
				var dN124 = 1.0 / 4.0 * (1.0 - ksi3);

				var aMatrix = new double[,]
					{
					{ -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, -N5 ,0.0, 0.0 , -N6, 0.0 ,0.0 , -N7, 0.0, 0.0, -N8, 0.0, 0.0, N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12, 0.0, 0.0 },
					{ 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, -N5 ,0.0, 0.0 , -N6, 0.0 ,0.0 , -N7, 0.0, 0.0, -N8, 0.0,0.0, N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12, 0.0 },
					{ 0.0, 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, -N5 ,0.0, 0.0 , -N6, 0.0 ,0.0 , -N7, 0.0, 0.0, -N8, 0.0, 0.0, N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12 }
					};

				var da1Matrix = new double[,]
					{
					{ -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, -dN51 ,0.0, 0.0 , -dN61, 0.0 ,0.0 , -dN71, 0.0, 0.0, -dN81, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, -dN51 ,0.0, 0.0 , -dN61, 0.0 ,0.0 , -dN71, 0.0, 0.0, -dN81, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, -dN51 ,0.0, 0.0 , -dN61, 0.0 ,0.0 , -dN71, 0.0, 0.0, -dN81, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
					};

				var da11Matrix = new double[,]
					{
					{ -dN111 ,0.0, 0.0 ,-dN211 ,0.0 ,0.0 ,-dN311, 0.0, 0.0, -dN411, 0.0, 0.0, -dN511 ,0.0, 0.0 , -dN611, 0.0 ,0.0 , -dN711, 0.0, 0.0, -dN811, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN111 ,0.0, 0.0 ,-dN211 ,0.0 ,0.0 ,-dN311, 0.0, 0.0, -dN411, 0.0, 0.0, -dN511 ,0.0, 0.0 , -dN611, 0.0 ,0.0 , -dN711, 0.0, 0.0, -dN811, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, -dN111 ,0.0, 0.0 ,-dN211 ,0.0 ,0.0 ,-dN311, 0.0, 0.0, -dN411, 0.0, 0.0, -dN511 ,0.0, 0.0 , -dN611, 0.0 ,0.0 , -dN711, 0.0, 0.0, -dN811, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
					};

				var da12Matrix = new double[,]
					{
					{ -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, -dN512 ,0.0, 0.0 , -dN612, 0.0 ,0.0 , -dN712, 0.0, 0.0, -dN812, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, -dN512 ,0.0, 0.0 , -dN612, 0.0 ,0.0 , -dN712, 0.0, 0.0, -dN812, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, -dN512 ,0.0, 0.0 , -dN612, 0.0 ,0.0 , -dN712, 0.0, 0.0, -dN812, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
					};

				var da2Matrix = new double[,]
					{
					{ -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, -dN52 ,0.0, 0.0 , -dN62, 0.0 ,0.0 , -dN72, 0.0, 0.0, -dN82, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, -dN52 ,0.0, 0.0 , -dN62, 0.0 ,0.0 , -dN72, 0.0, 0.0, -dN82, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, -dN52 ,0.0, 0.0 , -dN62, 0.0 ,0.0 , -dN72, 0.0, 0.0, -dN82, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
					};

				var da22Matrix = new double[,]
					{
					{ -dN122 ,0.0, 0.0 ,-dN222 ,0.0 ,0.0 ,-dN322, 0.0, 0.0, -dN422, 0.0, 0.0, -dN522 ,0.0, 0.0 , -dN622, 0.0 ,0.0 , -dN722, 0.0, 0.0, -dN822, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN122 ,0.0, 0.0 ,-dN222 ,0.0 ,0.0 ,-dN322, 0.0, 0.0, -dN422, 0.0, 0.0, -dN522 ,0.0, 0.0 , -dN622, 0.0 ,0.0 , -dN722, 0.0, 0.0, -dN822, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, -dN122 ,0.0, 0.0 ,-dN222 ,0.0 ,0.0 ,-dN322, 0.0, 0.0, -dN422, 0.0, 0.0, -dN522 ,0.0, 0.0 , -dN622, 0.0 ,0.0 , -dN722, 0.0, 0.0, -dN822, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
					};

				var da3Matrix = new double[,]
					{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN93, 0.0, 0.0, dN103, 0.0, 0.0, dN113, 0.0, 0.0, dN123, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, dN93, 0.0, 0.0, dN103, 0.0, 0.0, dN113, 0.0, 0.0, dN123, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, dN93, 0.0, 0.0, dN103, 0.0, 0.0, dN113, 0.0, 0.0, dN123 }
					};

				var da33Matrix = new double[,]
					{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
					};
				var da4Matrix = new double[,]
					{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN94, 0.0, 0.0, dN104, 0.0, 0.0, dN114, 0.0, 0.0, dN124, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, dN94, 0.0, 0.0, dN104, 0.0, 0.0, dN114, 0.0, 0.0, dN124, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, dN94, 0.0, 0.0, dN104, 0.0, 0.0, dN114, 0.0, 0.0, dN124 }
					};

				var da44Matrix = new double[,]
					{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
					};

				var T1 = new Tuple<double[,], double[,], double[,], double[,], double[,]>(da1Matrix, da11Matrix, da2Matrix, da22Matrix, da12Matrix);
				var T2 = new Tuple<double[,], double[,], double[,], double[,]>(da3Matrix, da33Matrix, da4Matrix, da44Matrix);
				return new Tuple<double[,], Tuple<double[,], double[,], double[,], double[,], double[,]>, Tuple<double[,], double[,], double[,], double[,]>>(aMatrix, T1, T2);
			}
			else
			{
				var N1 = 1.0 / 4.0 * (1d + ksi1) * (1d + ksi2) - 1.0 / 4.0 * (1d + ksi1) * (1 - ksi2 * ksi2) - 1.0 / 4.0 * (1d - ksi1 * ksi1) * (1d + ksi2);
				var N2 = 1.0 / 4.0 * (1d - ksi1) * (1d + ksi2) - 1.0 / 4.0 * (1d - ksi1 * ksi1) * (1 + ksi2) - 1.0 / 4.0 * (1d - ksi1) * (1d - ksi2 * ksi2);
				var N3 = 1.0 / 4.0 * (1d - ksi1) * (1d - ksi2) - 1.0 / 4.0 * (1d - ksi1 * ksi1) * (1 - ksi2) - 1.0 / 4.0 * (1d - ksi1) * (1d - ksi2 * ksi2);
				var N4 = 1.0 / 4.0 * (1d + ksi1) * (1d - ksi2) - 1.0 / 4.0 * (1d - ksi1 * ksi1) * (1 - ksi2) - 1.0 / 4.0 * (1d + ksi1) * (1d - ksi2 * ksi2);
				var N5 = 1.0 / 2.0 * (1d - ksi1 * ksi1) * (1d + ksi2);
				var N6 = 1.0 / 2.0 * (1d - ksi1) * (1d - ksi2 * ksi2);
				var N7 = 1.0 / 2.0 * (1d - ksi1 * ksi1) * (1d - ksi2);
				var N8 = 1.0 / 2.0 * (1d + ksi1) * (1d - ksi2 * ksi2);

				var dN11 = 1.0 / 4.0 * (2d * ksi1 + ksi2) * (1d + ksi2);
				var dN21 = 1.0 / 4.0 * (2d * ksi1 - ksi2) * (1d + ksi2);
				var dN31 = 1.0 / 4.0 * (2d * ksi1 + ksi2) * (1d - ksi2);
				var dN41 = 1.0 / 4.0 * (2d * ksi1 - ksi2) * (1d - ksi2);
				var dN51 = -ksi1 * (1d + ksi2);
				var dN61 = -1.0 / 2.0 * (1d - ksi2 * ksi2);
				var dN71 = -ksi1 * (1d - ksi2);
				var dN81 = 1.0 / 2.0 * (1d - ksi2 * ksi2);

				var dN111 = 1.0 / 2.0 * (1d + ksi2);
				var dN211 = 1.0 / 2.0 * (1d + ksi2);
				var dN311 = 1.0 / 2.0 * (1d - ksi2);
				var dN411 = 1.0 / 2.0 * (1d - ksi2);
				var dN511 = -(1d + ksi2);
				var dN611 = 0d;
				var dN711 = (ksi2 - 1d);
				var dN811 = 0d;

				var dN12 = 1.0 / 4.0 * (ksi1 + 2d * ksi2) * (1d + ksi1);
				var dN22 = 1.0 / 4.0 * (-ksi1 + 2d * ksi2) * (1d - ksi1);
				var dN32 = 1.0 / 4.0 * (ksi1 + 2d * ksi2) * (1d - ksi1);
				var dN42 = 1.0 / 4.0 * (-ksi1 + 2d * ksi2) * (1d + ksi1);
				var dN52 = 1.0 / 2.0 * (1d - Math.Pow(ksi1, 2));
				var dN62 = -(1d - ksi1) * ksi2;
				var dN72 = -1.0 / 2.0 * (1d - Math.Pow(ksi1, 2));
				var dN82 = -(1d + ksi1) * ksi2;

				var dN122 = 1.0 / 2.0 * (1d + ksi1);
				var dN222 = 1.0 / 2.0 * (1d - ksi1);
				var dN322 = 1.0 / 2.0 * (1d - ksi1);
				var dN422 = 1.0 / 2.0 * (1d + ksi1);
				var dN522 = 0d;
				var dN622 = ksi1 - 1d;
				var dN722 = 0d;
				var dN822 = -1d - ksi1;

				var dN112 = 1.0 / 4.0 * (2d * ksi1 + 2d * ksi2 + 1d);
				var dN212 = 1.0 / 4.0 * (2d * ksi1 - 2d * ksi2 - 1d);
				var dN312 = 1.0 / 4.0 * (-2d * ksi1 - 2d * ksi2 + 1d);
				var dN412 = 1.0 / 4.0 * (-2d * ksi1 + 2d * ksi2 - 1d);
				var dN512 = -ksi1;
				var dN612 = ksi2;
				var dN712 = ksi1;
				var dN812 = -ksi2;

				var N9 = 1.0 / 4.0 * (1d + ksi3) * (1d + ksi4) - 1.0 / 4.0 * (1d + ksi3) * (1 - ksi4 * ksi4) - 1.0 / 4.0 * (1d - ksi3 * ksi3) * (1d + ksi4);
				var N10 = 1.0 / 4.0 * (1d - ksi3) * (1d + ksi4) - 1.0 / 4.0 * (1d - ksi3 * ksi3) * (1 + ksi4) - 1.0 / 4.0 * (1d - ksi3) * (1d - ksi4 * ksi4);
				var N11 = 1.0 / 4.0 * (1d - ksi3) * (1d - ksi4) - 1.0 / 4.0 * (1d - ksi3 * ksi3) * (1 - ksi4) - 1.0 / 4.0 * (1d - ksi3) * (1d - ksi4 * ksi4);
				var N12 = 1.0 / 4.0 * (1d + ksi3) * (1d - ksi4) - 1.0 / 4.0 * (1d - ksi3 * ksi3) * (1 - ksi4) - 1.0 / 4.0 * (1d + ksi3) * (1d - ksi4 * ksi4);
				var N13 = 1.0 / 2.0 * (1d - ksi3 * ksi3) * (1d + ksi4);
				var N14 = 1.0 / 2.0 * (1d - ksi3) * (1d - ksi4 * ksi4);
				var N15 = 1.0 / 2.0 * (1d - ksi3 * ksi3) * (1d - ksi4);
				var N16 = 1.0 / 2.0 * (1d + ksi3) * (1d - ksi4 * ksi4);

				var dN93 = 1.0 / 4.0 * (2d * ksi3 + ksi4) * (1d + ksi4);
				var dN103 = 1.0 / 4.0 * (2d * ksi3 - ksi4) * (1d + ksi4);
				var dN113 = 1.0 / 4.0 * (2d * ksi3 + ksi4) * (1d - ksi4);
				var dN123 = 1.0 / 4.0 * (2d * ksi3 - ksi4) * (1d - ksi4);
				var dN133 = -ksi3 * (1d + ksi4);
				var dN143 = -1.0 / 2.0 * (1d - ksi4 * ksi4);
				var dN153 = -ksi3 * (1d - ksi4);
				var dN163 = 1.0 / 2.0 * (1d - ksi4 * ksi4);

				var dN933 = 1.0 / 2.0 * (1d + ksi4);
				var dN1033 = 1.0 / 2.0 * (1d + ksi4);
				var dN1133 = 1.0 / 2.0 * (1d - ksi4);
				var dN1233 = 1.0 / 2.0 * (1d - ksi4);
				var dN1333 = -(1d + ksi4);
				var dN1433 = 0d;
				var dN1533 = (ksi4 - 1d);
				var dN1633 = 0d;

				var dN94 = 1.0 / 4.0 * (ksi3 + 2d * ksi4) * (1d + ksi3);
				var dN104 = 1.0 / 4.0 * (-ksi3 + 2d * ksi4) * (1d - ksi3);
				var dN114 = 1.0 / 4.0 * (ksi3 + 2d * ksi4) * (1d - ksi3);
				var dN124 = 1.0 / 4.0 * (-ksi3 + 2d * ksi4) * (1d + ksi3);
				var dN134 = 1.0 / 2.0 * (1d - Math.Pow(ksi3, 2));
				var dN144 = -(1d - ksi3) * ksi4;
				var dN154 = -1.0 / 2.0 * (1d - Math.Pow(ksi3, 2));
				var dN164 = -(1d + ksi3) * ksi4;

				var dN944 = 1.0 / 2.0 * (1d + ksi3);
				var dN1044 = 1.0 / 2.0 * (1d - ksi3);
				var dN1144 = 1.0 / 2.0 * (1d - ksi3);
				var dN1244 = 1.0 / 2.0 * (1d + ksi3);
				var dN1344 = 0d;
				var dN1444 = ksi3 - 1d;
				var dN1544 = 0d;
				var dN1644 = -1d - ksi3;

				var aMatrix = new double[,]
					{
					{ -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, -N5 ,0.0, 0.0 , -N6, 0.0 ,0.0 , -N7, 0.0, 0.0, -N8, 0.0, 0.0, N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12, 0.0, 0.0, N13, 0.0, 0.0, N14, 0.0, 0.0, N15, 0.0, 0.0, N16, 0.0, 0.0 },
					{ 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, -N5 ,0.0, 0.0 , -N6, 0.0 ,0.0 , -N7, 0.0, 0.0, -N8, 0.0,0.0, N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12, 0.0, 0.0, N13, 0.0, 0.0, N14, 0.0, 0.0, N15, 0.0, 0.0, N16, 0.0 },
					{ 0.0, 0.0, -N1 ,0.0, 0.0 ,-N2 ,0.0 ,0.0 ,-N3, 0.0, 0.0, -N4, 0.0, 0.0, -N5 ,0.0, 0.0 , -N6, 0.0 ,0.0 , -N7, 0.0, 0.0, -N8, 0.0, 0.0, N9, 0.0, 0.0, N10, 0.0, 0.0, N11, 0.0, 0.0, N12, 0.0, 0.0, N13, 0.0, 0.0, N14, 0.0, 0.0, N15, 0.0, 0.0, N16 }
					};

				var da1Matrix = new double[,]
					{
					{ -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, -dN51 ,0.0, 0.0 , -dN61, 0.0 ,0.0 , -dN71, 0.0, 0.0, -dN81, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, -dN51 ,0.0, 0.0 , -dN61, 0.0 ,0.0 , -dN71, 0.0, 0.0, -dN81, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, -dN11 ,0.0, 0.0 ,-dN21 ,0.0 ,0.0 ,-dN31, 0.0, 0.0, -dN41, 0.0, 0.0, -dN51 ,0.0, 0.0 , -dN61, 0.0 ,0.0 , -dN71, 0.0, 0.0, -dN81, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
					};

				var da11Matrix = new double[,]
					{
					{ -dN111 ,0.0, 0.0 ,-dN211 ,0.0 ,0.0 ,-dN311, 0.0, 0.0, -dN411, 0.0, 0.0, -dN511 ,0.0, 0.0 , -dN611, 0.0 ,0.0 , -dN711, 0.0, 0.0, -dN811, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN111 ,0.0, 0.0 ,-dN211 ,0.0 ,0.0 ,-dN311, 0.0, 0.0, -dN411, 0.0, 0.0, -dN511 ,0.0, 0.0 , -dN611, 0.0 ,0.0 , -dN711, 0.0, 0.0, -dN811, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  },
					{ 0.0, 0.0, -dN111 ,0.0, 0.0 ,-dN211 ,0.0 ,0.0 ,-dN311, 0.0, 0.0, -dN411, 0.0, 0.0, -dN511 ,0.0, 0.0 , -dN611, 0.0 ,0.0 , -dN711, 0.0, 0.0, -dN811, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
					};

				var da12Matrix = new double[,]
					{
					{ -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, -dN512 ,0.0, 0.0 , -dN612, 0.0 ,0.0 , -dN712, 0.0, 0.0, -dN812, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, -dN512 ,0.0, 0.0 , -dN612, 0.0 ,0.0 , -dN712, 0.0, 0.0, -dN812, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  },
					{ 0.0, 0.0, -dN112 ,0.0, 0.0 ,-dN212 ,0.0 ,0.0 ,-dN312, 0.0, 0.0, -dN412, 0.0, 0.0, -dN512 ,0.0, 0.0 , -dN612, 0.0 ,0.0 , -dN712, 0.0, 0.0, -dN812, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
					};

				var da2Matrix = new double[,]
					{
					{ -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, -dN52 ,0.0, 0.0 , -dN62, 0.0 ,0.0 , -dN72, 0.0, 0.0, -dN82, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, -dN52 ,0.0, 0.0 , -dN62, 0.0 ,0.0 , -dN72, 0.0, 0.0, -dN82, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, -dN12 ,0.0, 0.0 ,-dN22 ,0.0 ,0.0 ,-dN32, 0.0, 0.0, -dN42, 0.0, 0.0, -dN52 ,0.0, 0.0 , -dN62, 0.0 ,0.0 , -dN72, 0.0, 0.0, -dN82, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
					};

				var da22Matrix = new double[,]
					{
					{ -dN122 ,0.0, 0.0 ,-dN222 ,0.0 ,0.0 ,-dN322, 0.0, 0.0, -dN422, 0.0, 0.0, -dN522 ,0.0, 0.0 , -dN622, 0.0 ,0.0 , -dN722, 0.0, 0.0, -dN822, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, -dN122 ,0.0, 0.0 ,-dN222 ,0.0 ,0.0 ,-dN322, 0.0, 0.0, -dN422, 0.0, 0.0, -dN522 ,0.0, 0.0 , -dN622, 0.0 ,0.0 , -dN722, 0.0, 0.0, -dN822, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
					{ 0.0, 0.0, -dN122 ,0.0, 0.0 ,-dN222 ,0.0 ,0.0 ,-dN322, 0.0, 0.0, -dN422, 0.0, 0.0, -dN522 ,0.0, 0.0 , -dN622, 0.0 ,0.0 , -dN722, 0.0, 0.0, -dN822, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
					};

				var da3Matrix = new double[,]
					{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN93, 0.0, 0.0, dN103, 0.0, 0.0, dN113, 0.0, 0.0, dN123, 0.0, 0.0, dN133, 0.0, 0.0, dN143, 0.0, 0.0, dN153, 0.0, 0.0, dN163, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, dN93, 0.0, 0.0, dN103, 0.0, 0.0, dN113, 0.0, 0.0, dN123, 0.0, 0.0, dN133, 0.0, 0.0, dN143, 0.0, 0.0, dN153, 0.0, 0.0, dN163, 0.0 },
					{0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, dN93, 0.0, 0.0, dN103, 0.0, 0.0, dN113, 0.0, 0.0, dN123, 0.0, 0.0, dN133, 0.0, 0.0, dN143, 0.0, 0.0, dN153, 0.0, 0.0, dN163 }
					};

				var da33Matrix = new double[,]
					{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN933, 0.0, 0.0, dN1033, 0.0, 0.0, dN1133, 0.0, 0.0, dN1233, 0.0, 0.0, dN1333, 0.0, 0.0, dN1433, 0.0, 0.0, dN1533, 0.0, 0.0, dN1633, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, dN933, 0.0, 0.0, dN1033, 0.0, 0.0, dN1133, 0.0, 0.0, dN1233, 0.0, 0.0, dN1333, 0.0, 0.0, dN1433, 0.0, 0.0, dN1533, 0.0, 0.0, dN1633, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, dN933, 0.0, 0.0, dN1033, 0.0, 0.0, dN1133, 0.0, 0.0, dN1233, 0.0, 0.0, dN1333, 0.0, 0.0, dN1433, 0.0, 0.0, dN1533, 0.0, 0.0, dN1633 }
					};

				var da4Matrix = new double[,]
					{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN94, 0.0, 0.0, dN104, 0.0, 0.0, dN114, 0.0, 0.0, dN124, 0.0, 0.0, dN134, 0.0, 0.0, dN144, 0.0, 0.0, dN154, 0.0, 0.0, dN164, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, dN94, 0.0, 0.0, dN104, 0.0, 0.0, dN114, 0.0, 0.0, dN124, 0.0, 0.0, dN134, 0.0, 0.0, dN144, 0.0, 0.0, dN154, 0.0, 0.0, dN164, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, dN94, 0.0, 0.0, dN104, 0.0, 0.0, dN114, 0.0, 0.0, dN124, 0.0, 0.0, dN134, 0.0, 0.0, dN144, 0.0, 0.0, dN154, 0.0, 0.0, dN164 }
					};

				var da44Matrix = new double[,]
					{
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dN944, 0.0, 0.0, dN1044, 0.0, 0.0, dN1144, 0.0, 0.0, dN1244, 0.0, 0.0, dN1344, 0.0, 0.0, dN1444, 0.0, 0.0, dN1544, 0.0, 0.0, dN1644, 0.0, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, dN944, 0.0, 0.0, dN1044, 0.0, 0.0, dN1144, 0.0, 0.0, dN1244, 0.0, 0.0, dN1344, 0.0, 0.0, dN1444, 0.0, 0.0, dN1544, 0.0, 0.0, dN1644, 0.0 },
					{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0 , 0.0, 0.0 ,0.0 , 0.0, 0.0, 0.0, 0.0, 0.0,0.0, dN944, 0.0, 0.0, dN1044, 0.0, 0.0, dN1144, 0.0, 0.0, dN1244, 0.0, 0.0, dN1344, 0.0, 0.0, dN1444, 0.0, 0.0, dN1544, 0.0, 0.0, dN1644 }
					};

				var T1 = new Tuple<double[,], double[,], double[,], double[,], double[,]>(da1Matrix, da11Matrix, da2Matrix, da22Matrix, da12Matrix);
				var T2 = new Tuple<double[,], double[,], double[,], double[,]>(da3Matrix, da33Matrix, da4Matrix, da44Matrix);
				return new Tuple<double[,], Tuple<double[,], double[,], double[,], double[,], double[,]>, Tuple<double[,], double[,], double[,], double[,]>>(aMatrix, T1, T2);
			}
		}

		private Tuple<Tuple<double[], double[], double[], double[], double[]>, double[,], double, double[,], double[], double[,], double[,]> MasterSurfaceGeometry(double[,] da1Matrix, double[,] da2Matrix, double[,] da11Matrix,
			double[,] da12Matrix, double[,] da22Matrix)
		{
			var xupd = NodalXUpdated().Scale(-1.0);
			var surfaceVector1 = Matrix.CreateFromArray(da1Matrix).Multiply(xupd);
			var surfaceVector2 = Matrix.CreateFromArray(da2Matrix).Multiply(xupd);
			var surfaceVectorDerivative11 = Matrix.CreateFromArray(da11Matrix).Multiply(xupd);
			var surfaceVectorDerivative12 = Matrix.CreateFromArray(da12Matrix).Multiply(xupd);
			var surfaceVectorDerivative22 = Matrix.CreateFromArray(da22Matrix).Multiply(xupd);
			var surafaceVectors = new Tuple<double[], double[], double[], double[], double[]>(surfaceVector1, surfaceVector2,
				surfaceVectorDerivative11, surfaceVectorDerivative12, surfaceVectorDerivative22);
			var m = new double[,]
					{
						{ surfaceVector1.DotProduct(surfaceVector1), surfaceVector1.DotProduct(surfaceVector2) },
						{ surfaceVector2.DotProduct(surfaceVector1), surfaceVector2.DotProduct(surfaceVector2) }
					};
			var detm = m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0];
			var mInv = new double[,]
					{
						{ m[1,1] / detm, - m[0,1] / detm },
						{ -m[1,0] / detm, m[0,0] / detm }
					};
			var normalUnitVec = (surfaceVector1.CrossProduct(surfaceVector2)).Scale(1.0 / Math.Pow(detm, 0.5));
			var curvatureTensor = new double[,]
					{
						{ 0.0, 0.0 },
						{ 0.0, 0.0 }
					};
			var curvatureContravariantTensor = new double[,]
					{
						{ 0.0, 0.0 },
						{ 0.0, 0.0 }
					};
			if (MasterSurfaceOrder >= 2)
			{
				curvatureTensor[0, 0] = surfaceVectorDerivative11.DotProduct(normalUnitVec);
				curvatureTensor[0, 1] = surfaceVectorDerivative12.DotProduct(normalUnitVec);
				curvatureTensor[1, 0] = surfaceVectorDerivative12.DotProduct(normalUnitVec);
				curvatureTensor[1, 1] = surfaceVectorDerivative22.DotProduct(normalUnitVec);
				curvatureContravariantTensor[0, 0] = curvatureTensor[0, 0] * mInv[0, 0] * mInv[0, 0] +
					curvatureTensor[0, 1] * mInv[0, 0] * mInv[1, 0] +
					curvatureTensor[1, 0] * mInv[0, 1] * mInv[0, 0] +
					curvatureTensor[1, 1] * mInv[0, 1] * mInv[1, 0];
				curvatureContravariantTensor[0, 1] = curvatureTensor[0, 0] * mInv[0, 0] * mInv[0, 1] +
					curvatureTensor[0, 1] * mInv[0, 0] * mInv[1, 1] +
					curvatureTensor[1, 0] * mInv[0, 1] * mInv[0, 1] +
					curvatureTensor[1, 1] * mInv[0, 1] * mInv[1, 1];
				curvatureContravariantTensor[1, 0] = curvatureTensor[0, 0] * mInv[1, 0] * mInv[0, 0] +
					curvatureTensor[0, 1] * mInv[1, 0] * mInv[1, 0] +
					curvatureTensor[1, 0] * mInv[1, 1] * mInv[0, 0] +
					curvatureTensor[1, 1] * mInv[1, 1] * mInv[1, 0];
				curvatureContravariantTensor[1, 1] = curvatureTensor[0, 0] * mInv[1, 0] * mInv[0, 1] +
					curvatureTensor[0, 1] * mInv[1, 0] * mInv[1, 1] +
					curvatureTensor[1, 0] * mInv[1, 1] * mInv[0, 1] +
					curvatureTensor[1, 1] * mInv[1, 1] * mInv[1, 1];
			}
			return new Tuple<Tuple<double[], double[], double[], double[], double[]>, double[,], double, double[,], double[], double[,], double[,]>(surafaceVectors, m, detm, mInv, normalUnitVec, curvatureTensor, curvatureContravariantTensor);
		}
		private double SlaveSurfaceMetricDet(double[,] da3Matrix, double[,] da4Matrix)
		{
			var xupd = NodalXUpdated();
			var surfaceVector1 = Matrix.CreateFromArray(da3Matrix).Multiply(xupd);
			var surfaceVector2 = Matrix.CreateFromArray(da4Matrix).Multiply(xupd);
			var m = new double[,]
			{
				{ surfaceVector1.DotProduct(surfaceVector1), surfaceVector1.DotProduct(surfaceVector2) },
				{ surfaceVector2.DotProduct(surfaceVector1), surfaceVector2.DotProduct(surfaceVector2) }
			};
			var detm = m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0];
			return detm;
		}
		private double CalculatePenetration(double[,] aMatrix, double[] n)
		{
			var AT = Matrix.CreateFromArray(aMatrix).Transpose();
			var AT_n = AT.Multiply(n);
			var xupd = NodalXUpdated();
			var normalGap = xupd.DotProduct(AT_n);
			return normalGap;
		}
		private double[] CalculateDeltaKsi(double[] masterSlaveRelativeVector, double[,] aMatrix, double[,] metricTensor,
			double[] surfaceVector1, double[] surfaceVector2, double[] surfaceVectorDerivative11, double[] surfaceVectorDerivative12,
			double[] surfaceVectorDerivative22, double detm, double[] xupd)
		{
			var deltaKsi = new double[2];
			if (MasterSurfaceOrder == 1)
			{
				var e = surfaceVectorDerivative12.DotProduct(Matrix.CreateFromArray(aMatrix).Multiply(xupd));
				var f = new double[]
				{
					surfaceVector1.DotProduct(Matrix.CreateFromArray(aMatrix).Multiply(xupd)),
					surfaceVector2.DotProduct(Matrix.CreateFromArray(aMatrix).Multiply(xupd))
				};
				var mat = new double[,]
				{
					{ metricTensor[1,1], e - metricTensor[0,1] },
					{ e - metricTensor[1,0], metricTensor[0,0] }
				};
				var scalar = 1.0 / (detm - Math.Pow(e, 2) + 2 * e * metricTensor[0, 1]);
				deltaKsi = (Matrix.CreateFromArray(mat).Multiply(f)).Scale(scalar);
			}
			else
			{
				var detDDF = (surfaceVector1.DotProduct(surfaceVector1) - surfaceVectorDerivative11.DotProduct(masterSlaveRelativeVector)) *
								(surfaceVector2.DotProduct(surfaceVector2) - surfaceVectorDerivative22.DotProduct(masterSlaveRelativeVector)) -
								(surfaceVector1.DotProduct(surfaceVector2) - surfaceVectorDerivative12.DotProduct(masterSlaveRelativeVector)) *
								(surfaceVector2.DotProduct(surfaceVector1) - surfaceVectorDerivative12.DotProduct(masterSlaveRelativeVector));
				var scalar = -1.0 / detDDF;
				var mat = new double[,]
				{
					{
						surfaceVector2.DotProduct(surfaceVector2) - surfaceVectorDerivative22.DotProduct(masterSlaveRelativeVector),
						surfaceVectorDerivative12.DotProduct(masterSlaveRelativeVector) - surfaceVector1.DotProduct(surfaceVector2)
					},
					{
						surfaceVectorDerivative12.DotProduct(masterSlaveRelativeVector) - surfaceVector2.DotProduct(surfaceVector1),
						surfaceVector1.DotProduct(surfaceVector1) - surfaceVectorDerivative11.DotProduct(masterSlaveRelativeVector)
					}
				};
				var vector = new double[]
				{
					surfaceVector1.DotProduct(masterSlaveRelativeVector),
					surfaceVector2.DotProduct(masterSlaveRelativeVector)
				};
				vector = vector.Scale(-1.0);
				deltaKsi = (Matrix.CreateFromArray(mat).Multiply(vector)).Scale(scalar);
			}
			return deltaKsi;
		}
		private double[] Project(double ksi1Initial, double ksi2Initial, double ksi3, double ksi4)
		{
			var maxIterations = 1000;
			var tol = Math.Pow(10.0, -6.0);
			var deltaKsi = new double[2];
			var ksi = new double[] { ksi1Initial, ksi2Initial };
			var xUpdated = NodalXUpdated();
			for (var i = 1; i <= maxIterations; i++)
			{
				var aMatrices = CalculatePositionMatrix(ksi[0], ksi[1], ksi3, ksi4);
				var masterSlaveRelativeVector = Matrix.CreateFromArray(aMatrices.Item1).Multiply(xUpdated);
				var masterSurfaceGeometry = MasterSurfaceGeometry(aMatrices.Item2.Item1, aMatrices.Item2.Item3, aMatrices.Item2.Item2,
																aMatrices.Item2.Item5, aMatrices.Item2.Item4);
				deltaKsi = CalculateDeltaKsi(masterSlaveRelativeVector, aMatrices.Item1, masterSurfaceGeometry.Item2,
					masterSurfaceGeometry.Item1.Item1, masterSurfaceGeometry.Item1.Item2, masterSurfaceGeometry.Item1.Item3,
					masterSurfaceGeometry.Item1.Item4, masterSurfaceGeometry.Item1.Item5, masterSurfaceGeometry.Item3,
					xUpdated);
				ksi[0] += deltaKsi[0];
				ksi[1] += deltaKsi[1];
				if (deltaKsi.Norm2() <= tol)
				{
					break;
				}
			}
			if (deltaKsi.Norm2() > tol)
			{
				throw new Exception("CPP not found in current iterations");
			}
			else
			{
				return ksi;
			}
		}

		private Tuple<double[], double[]> GaussPoints()
		{
			var iP = IntegrationPointsPerNaturalAxis;
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
			else if (iP == 7)
			{
				gaussPoints[0] = -0.94910;
				gaussPoints[1] = -0.74153;
				gaussPoints[2] = -0.40584;
				gaussPoints[3] = 0.0;
				gaussPoints[4] = 0.40584;
				gaussPoints[5] = 0.74153;
				gaussPoints[6] = 0.94910;
				gaussWeights[0] = 0.12948;
				gaussWeights[1] = 0.27970;
				gaussWeights[2] = 0.38183;
				gaussWeights[3] = 0.41795;
				gaussWeights[4] = 0.38183;
				gaussWeights[5] = 0.27970;
				gaussWeights[6] = 0.12948;
			}
			else if (iP == 8)
			{
				gaussPoints[0] = -0.96028;
				gaussPoints[1] = -0.79666;
				gaussPoints[2] = -0.52553;
				gaussPoints[3] = -0.18343;
				gaussPoints[4] = 0.18343;
				gaussPoints[5] = 0.52553;
				gaussPoints[6] = 0.79666;
				gaussPoints[7] = 0.96028;
				gaussWeights[0] = 0.10122;
				gaussWeights[1] = 0.22238;
				gaussWeights[2] = 0.31370;
				gaussWeights[3] = 0.36268;
				gaussWeights[4] = 0.36268;
				gaussWeights[5] = 0.31370;
				gaussWeights[6] = 0.22238;
				gaussWeights[7] = 0.10122;
			}
			else if (iP == 9)
			{
				gaussPoints[0] = -0.96816;
				gaussPoints[1] = -0.83603;
				gaussPoints[2] = -0.61337;
				gaussPoints[3] = -0.32425;
				gaussPoints[4] = 0.0;
				gaussPoints[5] = 0.32425;
				gaussPoints[6] = 0.61337;
				gaussPoints[7] = 0.83603;
				gaussPoints[8] = 0.96816;
				gaussWeights[0] = 0.08127;
				gaussWeights[1] = 0.18064;
				gaussWeights[2] = 0.26061;
				gaussWeights[3] = 0.31234;
				gaussWeights[4] = 0.33023;
				gaussWeights[5] = 0.31234;
				gaussWeights[6] = 0.26061;
				gaussWeights[7] = 0.18064;
				gaussWeights[8] = 0.08127;
			}
			else if (iP == 10)
			{
				gaussPoints[0] = -0.97390;
				gaussPoints[1] = -0.86506;
				gaussPoints[2] = -0.67940;
				gaussPoints[3] = -0.43339;
				gaussPoints[4] = -0.14887;
				gaussPoints[5] = 0.14887;
				gaussPoints[6] = 0.43339;
				gaussPoints[7] = 0.67940;
				gaussPoints[8] = 0.86506;
				gaussPoints[9] = 0.97390;
				gaussWeights[0] = 0.06667;
				gaussWeights[1] = 0.14945;
				gaussWeights[2] = 0.21908;
				gaussWeights[3] = 0.26926;
				gaussWeights[4] = 0.29552;
				gaussWeights[5] = 0.29552;
				gaussWeights[6] = 0.26926;
				gaussWeights[7] = 0.21908;
				gaussWeights[8] = 0.14945;
				gaussWeights[9] = 0.06667;
			}
			return new Tuple<double[], double[]>(gaussPoints, gaussWeights);
		}

		#region IElementType Members

		public int ID { get; set; }
		public int SubdomainID { get; set; }
		public IReadOnlyList<INode> Nodes { get; }

		public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofs;

		public IReadOnlyList<INode> GetNodesForMatrixAssembly() => Nodes;

		private Matrix CalculateMainStiffnessPart(Matrix A, double[] n)
		{
			var nxn = n.TensorProduct(n);
			var nxn_A = nxn * A;
			var AT_nxn_A = A.Transpose() * nxn_A;
			var mainStiffnessMatrix = AT_nxn_A.Scale(penaltyFactor);
			return mainStiffnessMatrix;
		}

		private Matrix CalculateRotationalStiffnessPart(Matrix A, Matrix dA1, Matrix dA2, double[] n, double ksi3, Matrix mInv, double[] surfaceVector1, double[] surfaceVector2)
		{
			var coef1 = penaltyFactor * ksi3 * mInv[0, 0];
			var n_x_surfaceVector1 = n.TensorProduct(surfaceVector1);
			var surfaceVector1_x_n = surfaceVector1.TensorProduct(n);
			var firstTerm = dA1.Transpose() * n_x_surfaceVector1 * A;
			var secondTerm = A.Transpose() * surfaceVector1_x_n * dA1;
			var rotationalPart1 = (firstTerm + secondTerm).Scale(coef1);
			var coef2 = penaltyFactor * ksi3 * mInv[1, 0];
			var n_x_surfaceVector2 = n.TensorProduct(surfaceVector2);
			var firstTerm2 = dA1.Transpose() * n_x_surfaceVector2 * A;
			var secondTerm2 = A.Transpose() * surfaceVector1_x_n * dA2;
			var rotationalPart2 = (firstTerm2 + secondTerm2).Scale(coef2);
			var coef3 = penaltyFactor * ksi3 * mInv[0, 1];
			var surfaceVector2_x_n = surfaceVector2.TensorProduct(n);
			var firstTerm3 = dA2.Transpose() * n_x_surfaceVector1 * A;
			var secondTerm3 = A.Transpose() * surfaceVector2_x_n * dA1;
			var rotationalPart3 = (firstTerm3 + secondTerm3).Scale(coef3);
			var coef4 = penaltyFactor * ksi3 * mInv[1, 1];
			var firstTerm4 = dA2.Transpose() * n_x_surfaceVector2 * A;
			var secondTerm4 = A.Transpose() * surfaceVector2_x_n * dA2;
			var rotationalPart4 = (firstTerm4 + secondTerm4).Scale(coef4);
			var rotationalPart = rotationalPart1 + rotationalPart2 + rotationalPart3 + rotationalPart4;
			return rotationalPart;
		}
		private Matrix CalculateCurvatureStiffnessPart(Matrix A, double ksi3, Matrix h, double[] surfaceVector1, double[] surfaceVector2)
		{
			var coef1 = penaltyFactor * ksi3 * h[0, 0];
			var coef2 = penaltyFactor * ksi3 * h[1, 0];
			var coef3 = penaltyFactor * ksi3 * h[0, 1];
			var coef4 = penaltyFactor * ksi3 * h[1, 1];

			var surfaceVector1_x_surfaceVector1 = surfaceVector1.TensorProduct(surfaceVector1);
			var surfaceVector1_x_surfaceVector2 = surfaceVector1.TensorProduct(surfaceVector2);
			var surfaceVector2_x_surfaceVector1 = surfaceVector2.TensorProduct(surfaceVector1);
			var surfaceVector2_x_surfaceVector2 = surfaceVector2.TensorProduct(surfaceVector2);
			var firstTerm = (A.Transpose() * surfaceVector1_x_surfaceVector1 * A).Scale(coef1);
			var secondTerm = (A.Transpose() * surfaceVector2_x_surfaceVector1 * A).Scale(coef2);
			var thirdTerm = (A.Transpose() * surfaceVector1_x_surfaceVector2 * A).Scale(coef3);
			var fourthTerm = (A.Transpose() * surfaceVector2_x_surfaceVector2 * A).Scale(coef4);
			var curvaturePart = firstTerm + secondTerm + thirdTerm + fourthTerm;
			return curvaturePart;
		}
		public IMatrix StiffnessMatrix()
		{
			var globalStifnessMatrix = Matrix.CreateFromArray(new double[3 * Nodes.Count, 3 * Nodes.Count]);
			var gPArray = GaussPoints().Item1;
			var gWArray = GaussPoints().Item2;
			for (var i = 0; i < IntegrationPointsPerNaturalAxis; i++)
			{
				var ihta1 = gPArray[i];
				var gW1 = gWArray[i];
				for (var j = 0; j < IntegrationPointsPerNaturalAxis; j++)
				{
					var ihta2 = gPArray[j];
					var gW2 = gWArray[j];
					var ksi = Project(0.0, 0.0, ihta1, ihta2);
					if (Math.Abs(ksi[0]) <= 1.05 && Math.Abs(ksi[1]) <= 1.05)
					{
						var positionMatrices = CalculatePositionMatrix(ksi[0], ksi[1], ihta1, ihta2);
						var aMatrix = positionMatrices.Item1;
						var da1Matrix = positionMatrices.Item2.Item1;
						var da2Matrix = positionMatrices.Item2.Item3;
						var da11Matrix = positionMatrices.Item2.Item2;
						var da12Matrix = positionMatrices.Item2.Item5;
						var da22Matrix = positionMatrices.Item2.Item4;
						var masterSurfaceCharacteristics = MasterSurfaceGeometry(da1Matrix, da2Matrix, da11Matrix,
																				da12Matrix, da22Matrix);
						var m = masterSurfaceCharacteristics.Item4;
						var dRho1 = masterSurfaceCharacteristics.Item1.Item1;
						var dRho2 = masterSurfaceCharacteristics.Item1.Item2;
						var n = masterSurfaceCharacteristics.Item5;
						//n = NormalVectorChecked(n);
						var h = masterSurfaceCharacteristics.Item7;
						var ksi3 = CalculatePenetration(aMatrix, n);
						if (ksi3 <= 0)
						{
							var slaveMetricTensorDet = SlaveSurfaceMetricDet(positionMatrices.Item3.Item1,
								positionMatrices.Item3.Item3);
							var mainPart = CalculateMainStiffnessPart(Matrix.CreateFromArray(aMatrix), n);
							var rotationalPart = CalculateRotationalStiffnessPart(Matrix.CreateFromArray(aMatrix), Matrix.CreateFromArray(da1Matrix),
								Matrix.CreateFromArray(da2Matrix), n, ksi3,
								Matrix.CreateFromArray(m), dRho1, dRho2);
							var curvaturePart = Matrix.CreateFromArray(new double[3 * Nodes.Count, 3 * Nodes.Count]);
							if (MasterSurfaceOrder != 1)
							{
								curvaturePart = CalculateCurvatureStiffnessPart(Matrix.CreateFromArray(aMatrix), ksi3,
								Matrix.CreateFromArray(h), dRho1, dRho2);
							}
							var scalar = Math.Pow(slaveMetricTensorDet, 0.5) * gW1 * gW2;
							var integrationPointStifnessContribution = (mainPart + rotationalPart + curvaturePart) * scalar;
							globalStifnessMatrix += integrationPointStifnessContribution;
						}
					}
				}
			}
			//if (activeElements.Contains(ID))
			//{
			//	(new MGroup.LinearAlgebra.Output.Array2DWriter()).WriteToFile(globalStifnessMatrix.CopytoArray2D(), $@"C:\Users\Public\Documents\MSolve_output\STSContactStiffnessElemID{ID}_TimeStep{AnalysisState.newmarkIncrementNumber}_Iteration_{AnalysisState.loadControlIteration}.txt");

			//}
			return dofEnumerator.GetTransformedMatrix(globalStifnessMatrix);
		}

		public IMatrix PhysicsMatrix() => StiffnessMatrix();

		public IMatrix MassMatrix()
		{
			double[,] massMatrix = new double[3 * Nodes.Count, 3 * Nodes.Count];
			return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(massMatrix));
		}

		public IMatrix DampingMatrix()
		{
			double[,] dampingMatrix = new double[3 * Nodes.Count, 3 * Nodes.Count];
			return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(dampingMatrix));
		}

		public double[] CreateInternalGlobalForcesVector()
		{
			var internalGlobalForcesVector = new double[3 * Nodes.Count];
			var gPArray = GaussPoints().Item1;
			var gWArray = GaussPoints().Item2;
			for (var i = 0; i < IntegrationPointsPerNaturalAxis; i++)
			{
				var ihta1 = gPArray[i];
				var gW1 = gWArray[i];
				for (var j = 0; j < IntegrationPointsPerNaturalAxis; j++)
				{
					var ihta2 = gPArray[j];
					var gW2 = gWArray[j];
					var ksi = Project(0.0, 0.0, ihta1, ihta2);
					if (Math.Abs(ksi[0]) <= 1.05 && Math.Abs(ksi[1]) <= 1.05)
					{
						var positionMatrices = CalculatePositionMatrix(ksi[0], ksi[1], ihta1, ihta2);
						var aMatrix = positionMatrices.Item1;
						var da1Matrix = positionMatrices.Item2.Item1;
						var da2Matrix = positionMatrices.Item2.Item3;
						var da11Matrix = positionMatrices.Item2.Item2;
						var da12Matrix = positionMatrices.Item2.Item5;
						var da22Matrix = positionMatrices.Item2.Item4;
						var masterSurfaceCharacteristics = MasterSurfaceGeometry(da1Matrix, da2Matrix, da11Matrix,
																					da12Matrix, da22Matrix);
						var n = masterSurfaceCharacteristics.Item5;
						var ksi3 = CalculatePenetration(aMatrix, n);
						if (ksi3 <= 0)
						{
							var slaveMetricTensor = SlaveSurfaceMetricDet(positionMatrices.Item3.Item1, positionMatrices.Item3.Item3);
							var scalar = Math.Pow(slaveMetricTensor, 0.5) * gW1 * gW2;
							var AT = Matrix.CreateFromArray(aMatrix).Transpose();
							var AT_n = AT.Multiply(n);
							var internalForcesVectorIntegarationPointContribution = AT_n.Scale(penaltyFactor * ksi3 * scalar);
							internalGlobalForcesVector = internalGlobalForcesVector.Add(internalForcesVectorIntegarationPointContribution);
						}
					}
				}
			}
			//if (activeElements.Contains(ID))
			//{
			//	(new MGroup.LinearAlgebra.Output.Array1DWriter()).WriteToFile(internalGlobalForcesVector, $@"C:\Users\Public\Documents\MSolve_output\STSContactInternalForcesElemID{ID}_TimeStep{AnalysisState.newmarkIncrementNumber}_Iteration_{AnalysisState.loadControlIteration}.txt");

			//}
			return internalGlobalForcesVector;
		}

		public Tuple<double[], double[]> CalculateResponse(double[] local_Displacements)
		{
			//if (activeElements.Contains(ID))
			//{
			//	(new MGroup.LinearAlgebra.Output.Array1DWriter()).WriteToFile(local_Displacements, $@"C:\Users\Public\Documents\MSolve_output\STSContactTotalDisplacementElemID{ID}_TimeStep{AnalysisState.newmarkIncrementNumber}_Iteration_{AnalysisState.loadControlIteration}.txt");

			//}

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

		public double[] CalculateResponseIntegralForLogging(double[] localDisplacements) => CalculateResponseIntegral();

		public double[] CalculateResponseIntegral() => CreateInternalGlobalForcesVector();

		public void SaveConstitutiveLawState(IHaveState externalState) { }

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
