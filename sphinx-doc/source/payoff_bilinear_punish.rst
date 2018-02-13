.. payoff_bilinear_punish:

.. default-domain:: mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Punishment values and BLPs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The first step is to construct punishment values over each partition element
of the state space. This will turn out to be amenable to *separable bilinear
programming* formulations.

Punishment payoffs for the large player :math:`G` are constructed as follows.

* For each :math:`a\in A`, define a correspondence :math:`\tilde{{\mathcal{W}}}:D\times A \rightrightarrows \mathbb{R}^\mathcal{Z}` by

    .. math::
        :nowrap:
        
        \begin{align*}
        \tilde{{\mathcal{W}}}(\lambda',a)&:=\{w\in {\mathcal{W}}(\lambda')\,|\, \forall j\in \mathcal{Z}, \forall a'\in \tilde{A}
        \\
        & \delta[p^j(a')-p^j(a_j)]\cdot w\leq (1-\delta) [\phi(a')-\phi(a_j)]\}.
        \end{align*}
    
  .. note ::
    Note :math:`\tilde{\mathcal{W}}` encodes two requirements:
    
    *   Continuation values :math:`w` must be *consistent* with the set :math:`\tilde{\mathcal{W}}(\lambda')` where :math:`w` itself enforces the action-continuation-state profile pair :math:`(a,\lambda') \equiv (a,\lambda P(a))`; and

    * individually, the action-and-promised-value pair :math:`(a_j, w_j)`
      are optimal (i.e. "incentive compatible").
        

        
* Next, construct government punishment vectors 
    
  * :math:`(\check{\pi}_k)_{k=1}^K`, and 
  
  * :math:`(\hat{\pi}_k)_{k=1}^K`
    
  by letting
    .. math::   
        \pi(\lambda):=\max_{b\in B} \min_{a,\lambda',w}&[(1-\delta)\lambda\cdot v(a,b)+\delta \lambda'\cdot w],\\
        \text{s.t.}\,\,& \lambda'=\lambda P(a),\\
        & -\lambda b^T \geq 0,\\
        & w\in \tilde{{\mathcal{W}}}(\lambda',a),
        :label: pival_k
    
  Then define: :math:`\check{\pi}_k:=\min_{\lambda\in Q_k}\pi(\lambda)` and :math:`\hat{\pi}_k:=\max_{\lambda\in Q_k}\pi(\lambda)`.

.. note ::
    * If :math:`\mathcal{W}(\lambda')` is defined as a convex polytope, then :math:`\tilde{\mathcal{W}}(\lambda',a)` is also a convex polytope.
    * Given :math:`(a,b) \in A \times B`, the minimization program :eq:`pival_k` is a mild nonlinear programming problem--i.e. a *separable bilinear programming formulation*--of the following generic form
      
    .. math ::
        \min_{\lambda,w} & u^{T} + \lambda^{T} Q w
        \\
        s.t. & \lambda \in \mathcal{F}_1,
        \\
        & w \in \mathcal{F}_2.     
        :label: BPseparable
        
    * Moreover, :math:`\mathcal{F}_1` and :math:`\mathcal{F}_2` are disjoint convex and bounded polytopes.

.. topic:: Existence of Optimum

    If :math:`F_1` and :math:`F_2` are bounded then there exists an optimal solution of :eq:`BPseparable`, :math:`(\lambda^{\ast}, w^{\ast})`, such that :math:`\lambda^{\ast} \in \mathcal{F}_1` and :math:`w \in \mathcal{F}_2`.

    See [HPT2000]_ or [HT1996]_.

.. topic:: Separable BLP is NP-complete

    A bilinear program can be solved in NP time.
    
    See [Man1995]_.
    
|

Assume that we are solving the minimization problem in :eq:`BPseparable`. We employ a well-known deterministic global optimization algorithm known as branch-and-bound (BNB). First, a *relaxation* of the bilinear program is solved. Typically, this is done by solving inexpensive LPs. [CM2009]_ The solution of the relaxed problem yields a lower bound on the globally optimal solution which is a convex lower envelope. The relaxation is then iteratively refined, by refining the domain (feasible sets) and successively eliminating dominated local optima. (This is also a common method in solving integer linear programs.) An upper bound estimate of the optima can be found by using local nonlinear solvers (e.g. SNOPT and IPOPT) over each branch. Thus we have successively improved branching partitions of the domain (i.e. *branching*) and lower- and upper-bounding estimates (i.e. *bounding*) of the :math:`\epsilon`-global optimum. 

.. note ::
    The BNB algorithm we use follows [McC1976]_ and is implemented through the
    Bilinear Matrix Inequality BNB interface (BMIBNB) available in Stefan
    Lofberg's `YALMIP <http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Solvers.BMIBNB>`_. To solve the local lower-bounding LPs, we use the `GNU GLPK <http://www.gnu.org/software/glpk/>`_ open-source optimizer and to solve the upper-bounding nonlinear programs, we use either `SNOPT <http://www.sbsi-sol-optimize.com/asp/sol_products_snopt_desc.htm>`_ or `IPOPT <https://projects.coin-or.org/Ipopt>`_.

===================================
Implementing punishment values
===================================

The following pseudocode implements the punishment calculations for the
outer-approximation scheme:

.. topic:: Pseudocode

    **Input**: :math:`H, c, PolyLcons`

    | **For each** :math:`(Q_k,a,b) \in D \times A \times B`:

    * Markov map :math:`P \gets P(a)`

    * Simplex :math:`T \gets P(Q_k)`

    * Get :math:`\mathbf{I}(a,k) := \{ k' \in \mathbf{K}: k' \mapsto Q_{k'}
      \in D, Q_{k'} \cap T \neq \emptyset` \}
        
        * Uses: *xTriIndex* from **Simplex_IntersectPmap**

        | **For each** :math:`k' \in \mathbf{I}(a,k)`:

            * Get linear inequality representations of :math:`T \cap Q_{k'}` as

            .. math::

                (M_{k'}, d_{k'}) \gets xPolyLcons(a,k,k')
                    
            * Get current payoff profile :math:`v(a,b)`

            * Solve separable BLP 

            .. math::

               \check{\pi}(Q_k,a,b,k') & 
               = \min_{\lambda \in Q_k,w} \lambda [ (1-\delta)v(a,b) + \delta P w ]                  
               \\
               s.t. &
               \\
               M_{k'} (\lambda P)^T & \leq d_{k'}
               \\        
               \lambda b^T & \leq 0
               \\
               H w & \leq c^{o}(Q_k', \cdot)
               \\
               \delta[p^j(a')-p^j(a_j)]\cdot w & \leq (1-\delta) [\phi(a')-\phi(a_j)], \  \forall j \in \mathcal{Z}
            
    * Get :math:`\check{\pi}(Q_k,a,b) = \min_{k' \in \mathbf{I}(a,k)} \check{\pi}(Q_k,a,b,k')`.

 * Get :math:`\check{\pi}(Q_k)  = \max_{b \in B} \min_{a \in A}\check{\pi}(Q_k,a,b)`.


|

.. note ::
    
    * The *separable* constraint sets :math:`\mathcal{F}_1` and
      :math:`\mathcal{F}_2` for :math:`\lambda` and :math:`w`, respectively,
      are given by constraints in the BLP. These constraint say the following.
      
        * :math:`\mathcal{F}_1`: the first two constraints require :math:`\lambda \in Q_k` to be such that
            
            * for each :math:`a \in A, k \in \mathbf{K}` and :math:`k' \in
              \mathbf{I}(a,k)`, the resulting continuation state :math:`\lambda P(a) \in Q_{k'} \cap P(Q_k)`; and
            * given a fixed policy :math:`b`, the choice over :math:`\lambda`
              renders :math:`b` feasible according the the government budget
              constraint;
        * :math:`\mathcal{F}_2` is given by the requirements that :math:`w` be
            
            * consistent with respect to the step correspondence slice
              :math:`\mathcal{W}(k')` which has constant levels over the partition
              element :math:`Q_{k'}`; and
            * such that :math:`w` is incentive compatible for all small
              agents.
      
    * Constructing the punishment values :math:`\hat{\pi}_k` for the *inner-approximation* scheme is similar to what we did above in detail for :math:`\check{\pi}_k`. The only differences are
      
     #. in the second last step of the pseudocode above, replace that line with:
     
            .. math::        
                \hat{\pi}(Q_k,a,b) = \max_{k' \in \mathbf{I}(a,k)} \hat{\pi}(Q_k,a,b,k').             
        
        Of course we should also re-label all the :math:`\check{\pi}` notation for the punishment value function with :math:`\hat{\pi}`; and
        
     #. *maximize* over :math:`\lambda \in Q_k` in the main BLP problem.

|

.. topic:: Relevant functions :math:`\blacktriangleright`

    .. function:: Punish_Outer(self)

        Returns:

        * *pival* :
    
            * A :math:`(K \times 1)` numeric array containing elements
              :math:`\check{\pi}(Q_k)` where :math:`k \in \mathbf{K}`.

    .. seealso:: **PunishK**

    .. function:: Punish_Inner(self)

        Returns:

        * *pival* :
    
            * A :math:`(K \times 1)` numeric array containing elements
              :math:`\hat{\pi}(Q_k)` where :math:`k \in \mathbf{K}`.
        
    .. seealso:: **PunishK**



|
|
|

.. [CM2009] Carpara, Alberto and Michele Monaci. "Bidimensional packing by bilinear programming". *Mathematical Programming Series A*, 118, 75–108.

.. [HPT2000] Horst R, P. Pardalos and N. Thoai (2000): *Introduction to global optimization*. 2nd Edition, Boston: Springer, 2000.

.. [HT1996] Horst R, Hoang Tuy (1996): *Global Optimization*. 3rd Edition, New York: Springer.

.. [Man1995] Mangasarian, Olvi L. (1995): "The linear complementarity problem as a separable bilinear program". *Journal of Global Optimization*, 12, 1–7.


