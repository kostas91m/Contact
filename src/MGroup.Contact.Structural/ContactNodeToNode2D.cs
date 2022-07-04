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
	public class ContactNodeToNode2D : IStructuralElementType
	{
		private static readonly IDofType[] nodalDOFTypes = new IDofType[] { StructuralDof.TranslationX, StructuralDof.TranslationY };
		private static readonly IDofType[][] dofs = new IDofType[][] { nodalDOFTypes, nodalDOFTypes };
		private readonly double penaltyFactor;
		private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();
		private double[] DisplacementVector { get; set; }
		private double ContactArea { get; }

		public ContactNodeToNode2D(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplier,
			double contactArea)
		{
			this.penaltyFactor = penaltyFactorMultiplier * youngModulus;
			this.DisplacementVector = new double[4];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
		}
		public ContactNodeToNode2D(IReadOnlyList<INode> nodes, double penaltyFactor, double contactArea)
		{
			this.penaltyFactor = penaltyFactor;
			this.DisplacementVector = new double[4];
			this.ContactArea = contactArea;
			this.Nodes = nodes;
		}
		public ContactNodeToNode2D(IReadOnlyList<INode> nodes, double youngModulus, double penaltyFactorMultiplier,
			IElementDofEnumerator dofEnumerator)
			: this(nodes, youngModulus, penaltyFactorMultiplier, 1d)
		{
			this.dofEnumerator = dofEnumerator;
		}
		public ContactNodeToNode2D(IReadOnlyList<INode> nodes, double penaltyFactor, IElementDofEnumerator dofEnumerator)
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

		private double CalculateNormalGap()
		{
			double[,] A = CalculatePositionMatrix();
			var AT = Matrix.CreateFromArray(A).Transpose();
			double[] n = CalculateNormalUnitVector();
			var AT_n = Vector.CreateFromArray(AT.Multiply(n));
			var xupd = Vector.CreateFromArray(new double[]
			{
				Nodes[0].X + DisplacementVector[0],
				Nodes[0].Y + DisplacementVector[1],
				Nodes[1].X + DisplacementVector[2],
				Nodes[1].Y + DisplacementVector[3]
			});
			double normalGap = xupd.DotProduct(AT_n);
			return normalGap;
		}

		private double[] CalculateNormalUnitVector()
		{
			double X1 = Nodes[0].X;
			double Y1 = Nodes[0].Y;
			double X2 = Nodes[1].X;
			double Y2 = Nodes[1].Y;
			var normalVector = Vector.CreateFromArray(new double[] { X2 - X1, Y2 - Y1 });
			double normalVectorLength = normalVector.Norm2();
			double[] normalUnitVec = new double[] { normalVector[0] / normalVectorLength, normalVector[1] / normalVectorLength };
			return normalUnitVec;
		}

		private double[,] CalculatePositionMatrix()
		{
			double[,] aMatrix = new double[,]
				{
					{ -1,0,1,0},
					{0,-1,0,1 }
				};
			return aMatrix;
		}

		#region IElementType Members

		public int ID { get; set; }
		public int SubdomainID { get; set; }
		public IReadOnlyList<INode> Nodes { get; }

		public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

		public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes() => dofs;

		public IReadOnlyList<INode> GetNodesForMatrixAssembly() => Nodes;

		public IMatrix StiffnessMatrix()
		{
			double penetration = CalculateNormalGap();
			if (penetration <= 0)
			{
				var n = Vector.CreateFromArray(CalculateNormalUnitVector());
				var A = Matrix.CreateFromArray(CalculatePositionMatrix());
				var AT = A.Transpose();
				var nxn = n.TensorProduct(n);
				var nxn_A = nxn.MultiplyRight(A);
				var AT_nxn_A = AT.MultiplyRight(nxn_A);
				var globalStiffnessMatrix = AT_nxn_A.Scale(penaltyFactor);
				return dofEnumerator.GetTransformedMatrix(globalStiffnessMatrix);
			}
			else
			{
				double[,] globalStifnessMatrix = new double[4, 4];
				return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(globalStifnessMatrix));
			}
		}
		public IMatrix PhysicsMatrix()
		{
			return StiffnessMatrix();
		}
		public IMatrix MassMatrix()
		{
			double[,] massMatrix = new double[4, 4];
			return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(massMatrix));
		}

		public IMatrix DampingMatrix() 
		{
			double[,] dampingMatrix = new double[4, 4];
			return dofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromArray(dampingMatrix));
		}

		private double[] CalculateResponseIntegral(double[] localdDisplacements)
		{
			if (localdDisplacements != null)
			{
				IMatrix stiffness = StiffnessMatrix();
				return stiffness.Multiply(localdDisplacements);
			}
			else
			{
				return null;
			}
		}

		public Tuple<double[], double[]> CalculateResponse(double[] local_Displacements)
		{
			// WARNING: 1) No strains are computed 2) localdDisplacements are not used.
			double[] strains = null;
			double[] forces = CalculateResponseIntegral(local_Displacements);
			double[] stresses = Array.ConvertAll(forces, x => x / ContactArea);
			if (DisplacementVector == null || DisplacementVector.Length != local_Displacements.Length)
			{
				DisplacementVector = new double[local_Displacements.Length];
			}

			Array.Copy(local_Displacements, DisplacementVector, local_Displacements.Length);

			return new Tuple<double[], double[]>(strains, stresses);
		}

		public double[] CalculateResponseIntegralForLogging(double[] localDisplacements)
		{
			return CalculateResponseIntegral();
		}

		public double[] CalculateResponseIntegral() => CalculateResponseIntegral(DisplacementVector);

		//public double[] CalculateAccelerationResponseIntegral(IElementType element, IList<MassAccelerationLoad> loads)
		//{
		//	var accelerations = new double[4];
		//	IMatrix massMatrix = MassMatrix(element);

		//	int index = 0;
		//	foreach (MassAccelerationLoad load in loads)
		//		foreach (IDofType[] nodalDOFTypes in dofs)
		//			foreach (IDofType dofType in nodalDOFTypes)
		//			{
		//				if (dofType == load.DOF) accelerations[index] += load.Amount;
		//				index++;
		//			}

		//	return massMatrix.Multiply(accelerations);
		//}

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
