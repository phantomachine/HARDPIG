.. compute_statespace:

.. default-domain:: mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
State Space Computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Now we describe the implementation for the key tasks involved so far. We will need to:

* compute state state partitions: :math:`D := \cup_{k=1}^{K}Q_k`;

* construct transitions from each subspace :math:`Q_k` into corresponding sets
  :math:`P(a)(Q_k) \subset D`, for every possible action profile :math:`a \in A`;

* record intersections :math:`P(a)(Q_k) \cap Q_k'`, for every :math:`k,k' \in \mathbf{K}`; and

* sample uniform points from each :math:`Q_k`, and check for nonempty samples that end up transitioning to each respective intersecting set :math:`P(a)(Q_k) \cap Q_k'`.

-----------------------
Storage
-----------------------

We only need to store:

#. each index :math:`k' \in \mathbf{I}(a,k) \subseteq \mathbf{K}`, which refers to some partition element(s) :math:`Q_{k'}` whose subset :math:`Poly_{k'}` is accessible from :math:`Q_{k}` given operator :math:`P(a)`.

    * This suffices for indexing the correct slices of equilibrium payoff sets over the corresponding subset :math:`Poly_{k'}` of the state space :math:`D`.

#. the finite number of vertices of each :math:`Poly_{k'}` and the corresponding linear (weak) inequality representation of each :math:`Poly_{k'}`. 

    * This will become apparent later when we solve *separable* bilinear programming problems where it involves optimizing over these subsets :math:`Poly_{k'}` (when constructing max-min punishment values in the game).

#. the sub-samples from all the Hit-and-Run realizations that belong to every partition element
   :math:`Q_k` which will end up in particular intersections summarized by
   each polytope :math:`Poly_{k'}`. 

-----------------------
Implementation
-----------------------

Since we have finite partition elements :math:`Q_{k}` and finite action
profile set :math:`A \ni a`, then we can enumerate and store all intersections
previously denoted by :math:`\{ Poly_{k'(k,a)}: \forall a \in A, k \in
\mathbf{K} \}` or equivalently by their index sets :math:`\{ \mathbf{I}(a,k):
(a,k) \in A \times \mathbf{K} \}`.


.. topic:: Pseudocode 

    | **For each** :math:`(a,k)\in A \times \mathbf{K}`:  

    * Get vertex representation of :math:`Q_k \in D`
    
    * Set :math:`P(a)` as :math:`P`
    
    * Get vertex representation :math:`T` from :math:`P(Q_k)`
    
    *  Simulate Hit-and-Run uniform realizations in simplex :math:`T`. Get 
        .. math ::
            X := \{ \lambda_n \}_{n=1}^{N_{sim}} \leftarrow RandomPolyFill(N_{sim},T)

    * Set :math:`i \leftarrow 0`

        | **For each** :math:`k' \in \mathbf{K}`:

        * Get vertex representation of :math:`Q_k' \in D`
        
        * Get intersection :math:`InterPoly = Q_{k'} \cap T` (a polytope) as:
            
        .. math ::
            InterPoly \leftarrow PolyBool(Q_{k'},T)

        * If :math:`InterPoly \neq \emptyset`:
            
            * Set :math:`i \leftarrow i + 1`
            
            * Store index to partition elements :math:`Q_{k'}` when :math:`InterPoly` is nonempty. Set 
           
            .. math ::
                TriIndex(a,k)(i) \gets k'

            * Store vertex data of polytope. Set: 
            
            .. math ::
                PolyVerts(a,k,k') \gets Interpoly

            * Store linear inequality representation of the same polytope. Set: 
            
            .. math ::
                PolyLcons(a,k,k') \gets Vert2Lcon(InterPoly)

            * Map Monte Carlo realizations :math:`X` under operator :math:`P`. Set:
            
            .. math ::
                Y \gets P(X)

            * List members of :math:`Y` that end up in :math:`InterPoly`. Set: 
            
            .. math ::
                :nowrap:

                \begin{eqnarray*}
                Y_{in} & \gets & InPolytope(Y, InterPoly)
                \\
                \{ PolyRand(a,k,k'), in \} & \gets & Y_{in}
                \end{eqnarray*}
      
            * Record all vectors :math:`\{\lambda_{n}\} \subseteq X` that induce :math:`Y_{in}` under map :math:`P(a)`:
            
            .. math ::
                \{PolyQjRand(a,k,k'), Index\} \gets X(in)
            
            * Store corresponding indices :math:`\{n : \lambda_n P(a)\in Y_{in}\}`:
            
            .. math ::
                PolyIndexQjRand(a,k,k') \gets Index 

    * **Return**: :math:`TriIndex, PolyVerts, PolyLcons, PolyRand, PolyQjRand, PolyIndexQjRand` 
 

.. note::
    Computationally, we only need to construct sets :math:`\{
    Poly_{k'(k,a)}: \forall a \in A, k \in \mathbf{K} \}` (e.g. *PolyVerts* and
    *PolyLcons* in the pseudocode above) or :math:`\{ \mathbf{I}(a,k): (a,k) \in A
    \times \mathbf{K} \}` (i.e. *TriIndex* in the pseudocode above) once beforehand.


.. topic:: Relevant functions :math:`\blacktriangleright`

    .. function:: Simplex_IntersectPmap(self)

        Returns 4 possible output:

        * *QP* :
    
            * a structure variable containing all others below.

        * *TriIndex* :
    
            * a cell array containing indices :math:`k'(a,k) \in \mathbf{I}(a,k) \subseteq \mathbf{K}` of partition elements that have non-empty intersections with each simplicial image :math:`P(a)(Q_{k})`.

        * *PolyVerts* :
    
            * a cell array, where each cell is an array with rows corresponding to vertices of :math:`Poly_{k'(a,k)}`, a polytope contained in the partition element :math:`Q_{k}`. Each cell element is consistent with the index :math:`k'(a,j) \in \mathbf{I}(a,k) \subseteq \mathbf{K}` stored in *TriIndex*.

        * *PolyLcons* :
    
            * is a set of linear (weak) inequality constraint representation of *PolyVerts*.

        * *PolyRand* :
          
            * Realizations of random vectors :math:`\lambda_{s}P(a) \in Q_{k'}` where :math:`k' \in \mathbf{I}(a,k)` and :math:`\lambda_{s} \sim \textbf{U}[Q_{k}]`---i.e. is uniformly drawn  from :math:`Q_{k}` according to a Hit-and-Run algorithm [ HYPERLINK TO ALGORITHM DESCRIPTION ], classified according to each *PolyVerts{a}{k}{k'}*.
              
        * *PolyQjRand* :

            * Inverse of *PolyRand*. Each *PolyQjRand{a}{k}{k'}* gives the set of :math:`\lambda_{s} \sim \textbf{U}[Q_{k}]`, where under action profile :math:`a`, the induced vector is :math:`P(a)(\lambda_{s}) \in Q_{k'}` and :math:`k' \in \mathbf{I}(a,k)`.


        * *PolyIndexQjRand* :

            * Each *PolyIndexQjRand{a}{j}* gives the set of indices :math:`\{
              s \in \{1,...,N_{sim}\} : s \mapsto \lambda_{s} \in Q_{k}, \lambda_{s}P(a) \in Q_{k'}  \text{ and } k' \in \mathbf{I}(a,k)  \}`. The number of *Monte Carlo simulations* of these uniform vectors subject to each convex body :math:`Q_k` has to be prespecified as :Math:`N_{sim}`.

    .. seealso:: **PolyBool**, **Simplex_Intersect**, **Vert2Lcon**, **RandPolyFill**, **cprnd**



