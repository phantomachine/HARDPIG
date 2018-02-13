.. payoff_compute:

.. default-domain:: mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Computing approximate SSE payoffs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Now we are ready to describe the computation of the approximate SSE payoff correspondence. The basic idea is from [JYC2003]_ who use *linear program* (LP) formulations as the approximation. Our extension illustrates that when we have probability distributions (with finite support) as state variables, the approximate SSE payoff correspondence can be constructed via *bilinear program* (BLP) formulations.

Notation:

* Let

  .. math::
    v_j(a,b) := u(c^b(j))-\phi(a_j)
* Given:
  
  * a vector of agent actions :math:`a`, 
    
  * a government policy vector :math:`b`, and

  * a vector of continuation payoffs :math:`w`, 
    
  the vector of agents' expected payoffs is defined by
    
    .. math::
        E(a,b,w):=((1-\delta)v_j(a,b)+\delta P^j(a_j)\cdot w)_{j\in \mathcal{Z}}.




.. _outer_concept:

====================================
Outer Approximation: Conceptual
====================================

We can now define the outer approximation :math:`\mathbf{B}^o(\mathcal{W})`. 

* For each search subgradient :math:`h_l\in H` and each partition element :math:`Q_k`, let

    .. math::
        c^{o+}_l(k):=\max_{(a,b) \in A \times B,\lambda\in Q_k, w}& [h_l \cdot E(a,b,w)],
        \\
        \text{s.t.}\,\,&\lambda'=\lambda P(a),
        \\
        & \lambda b^T \leq 0,
        \\
        & w\in \tilde{\mathcal{W}}(\lambda', a),
        \\
        & (1-\delta)\lambda \cdot v(a,b)+\delta \lambda' \cdot w\geq \check{\pi}_k,
        :label: OuterBLP


* Then define
    .. math::
        
        \bar{\omega}^{o+}_k (\lambda) := 
        \begin{cases}
        \bigcap_{l=1}^L \{z\,|\, h_l\cdot z\leq  c^{o+}_l(k)\},& \text{if} \ \lambda \in Q_k,\\
        \emptyset, & \text{otherwise}.
        \end{cases}
        
.. note ::
    
    Since :math:`A \times B` is a finite set of action profiles, we can evaluate the program :eq:`OuterBLP` as a special class of a nonlinear optimization problem--a *nonseparable* bilinear program (BLP)--for each fixed :math:`(a,b) \in A \times B`. Then we can maximize over the set :math:`A \times B`, by table look-up.

====================================
Outer Approximation: Implementation
====================================


Now we deal with implementing the idea in :ref:`outer_concept`. The
outer-approximation scheme to construct :math:`\mathbf{B}^o(\mathcal{W}^o)` in
the set of problems in :eq:`OuterBLP` is computable by following the
pseudocode below:

.. topic:: Pseudocode

    **Input**: :math:`H, c, \hat{\pi}, Poly`

    | **For each** :math:`(Q_k,h_l,a) \in D \times H \times A`:

    * Markov map :math:`P \gets P(a)`

    * Simplex :math:`T \gets P(Q_k)`

    * Get Hit-and-Run uniform draws constrained to be in :math:`Q_k`: :math:`X:= \{ \lambda_n \}_{n=1}^{N_{sim}} \gets RandomPolyFill(N_{sim}, T)`

    * Get feasible set :math:`F \gets F(B,Q_k) := \{ ( \lambda, b) \in Q_k \times B  : \lambda \in X \subseteq Q_k,  -\lambda b^T \geq \mathbf{0} \}`

    * Get :math:`\mathbf{I}(a,k) := \{ k' \in \mathbf{K}: k' \mapsto Q_{k'}
      \in D, Q_{k'} \cap T \neq \emptyset` \}
        
        * Uses: *TriIndex* from **Simplex_IntersectPmap**

        | **For each** :math:`k' \in \mathbf{I}(a,k)`:
        
        * Get relevant feasible policy set :math:`F(k'; a,k) \gets \{ (\lambda,b) \in F : \lambda P \in Q_{k'} \}`
        
            | **For each** :math:`(\lambda,b) \in F(k';a,k)`:
            
            * Get current payoff profile :math:`v(a,b)`

            * Solve conditional LP:

            .. math::

               c_{+}^{o}(Q_k,h_l,a,k',\lambda,b) & 
               = \max_{w} h_l [ (1-\delta)v(a,b) + \delta P w ]                  
               \\
               s.t. &
              \\             
               H w & \leq c^{o}(Q_{k'}, h_l)
               \\
               \delta[p^j(a')-p^j(a_j)]\cdot w & \leq (1-\delta) [\phi(a')-\phi(a_j)], \  \forall j \in \mathcal{Z}
               \\
               -\lambda \delta P w & 
               \leq \lambda [ (1-\delta)v(a,b) ] - \check{\pi}(Q_{k})
        
        * Get :math:`c_{+}^{o}(Q_k,h_l,a,k') = \max_{(\lambda,b) \in F(k';a,k)} c_{+}^{o}(Q_k,h_l,a, k',\lambda,b)`.
            
    * Get :math:`c_{+}^{o}(Q_k,h_l,a) = \max_{k' \in \mathbf{I}(a,k)} c_{+}^{o}(Q_k,h_l,a,k')`.

 * Get :math:`c_{+}^{o}(Q_k,h_l) = \max_{a \in A} c_{+}^{o}(Q_k,h_l,a)`.

|
|

.. note ::
    * In the pseudocode, we can see that for every fixed :math:`(Q_k,a,k') \in D
      \times A \times \mathbf{I}(a,k)` and every feasible :math:`b`, the nested family of programming problems are *nonseparable bilinear programs* (BLP) in the variables :math:`(\lambda, w)`.
    * The inner most loop thus implements our Monte Carlo approach to approximately
      solve for an :math:`\epsilon`-global solution to the *nonseparable* BLPs.
    * Conditional on each draw of :math:`\lambda`, this becomes a standard
      linear program (LP) in :math:`w` within each innermost loop of the pseudocode.
    * Given the set of subgradients :math:`H`, an outer-approximation update on the initial step
      correspondence :math:`\mathcal{W}^o`, is now sufficiently summarized by
      :math:`\mathcal{W}^{o}_{+} = \mathbf{B}^o(\mathcal{W}^o) \gets (H,c_{+}^{o}(Q_k,h_l),\mathcal{W}^o)`. 


|

.. topic:: Relevant functions :math:`\blacktriangleright`

    .. function:: Admit_Outer_LPset(self)

        Returns:

        * *Cnew* :
    
            * A :math:`(L \times K)` numeric array containing elements
              :math:`c_{+}^{o}(Q_k,h_l)` where :math:`k \in \mathbf{K}` and :math:`h_l \in
              H` are, respectively, a partition element of the correspondence
              domain :math:`D`, and, search subgradient in direction indexed
              by :math:`l \in \{1,...,|H|\}`.

    .. seealso:: **Punish_Outer**

====================================
Inner Approximation: Conceptual
====================================

We now define the *inner approximation* of the SSE value correspondence operator as :math:`\mathbf{B}^i(\mathcal{W})` below. 

* Denote :math:`\tilde{H}(k)` as a finite set of :math:`\tilde{L}(k)` spherical codes (to be used as *approximation subgradients*, where each
  element is :math:`\tilde{h}_l (k)` and :math:`\| \tilde{h}_l (k) \|_2 = 1` for all
  :math:`l = 1,...,\tilde{L}` and :math:`k \in \mathbf{K}`.

* Assume an initial inner step-correspondence approximation of some
  convex-valued and compact-graph correspondence  
  
.. math::
    \mathcal{W}^{i} := \bigcup_{Q_k \in D}  \bigcap_{\tilde{h}_l (k) \in \tilde{H}(k)} 
    \left\{ z \in \mathbb{R}^{\mathcal{Z}} :
    \tilde{h}_l (k) z \leq c(Q_k,\tilde{h}_l) \right\}.
  
* Define another finite set of fixed :math:`L` *search subgradients*, made up also of
  spherical codes, :math:`H`, just as in the outer approximation method above. [#footnote_subgradients_inner]_
 
* For each *search subgradient* :math:`h_l\in H` and each partition element
  :math:`Q_k`, let

    .. math::
        V^{i+}_l(k):=\min_{\lambda\in Q_k} \max_{(a,b) \in A \times B, w}& [h_l \cdot E(a,b,w)],
        \\
        \text{s.t.}\,\,&\lambda'=\lambda P(a),
        \\
        & \lambda b^T \leq 0,
        \\
        & w\in \tilde{\mathcal{W}^i}(\lambda',a),
        \\
        & (1-\delta)\lambda \cdot v(a,b)+\delta \lambda' \cdot w\geq \hat{\pi}_k,
        :label: InnerBLP

  Set :math:`V^{i+}_l(k) = -\infty` if the optimizer set is empty.
    
* In contrast to :ref:`outer_concept`, obtain the following additional
  step. 
    
    * Let :math:`(a_{l}^{\ast}(k), b_{l}^{\ast}(k), w_{l}^{\ast}(k))` denote the maximizers in direction :math:`h_l` and over domain partition element :math:`Q_k`, that induce the level :math:`V^{i+}_l(k)` above.

    * Then the corresponding vector of agent payoffs is
      
    .. math::
        z^{+}_l (k) := E(a_{l}^{\ast}(k), b_{l}^{\ast}(k), w_{l}^{\ast}(k)).

    * Define the set of vertices :math:`Z(k) = \{ z^{+}_l (k) : l =
      1,...,L \}` and let :math:`\mathcal{W}^{i+}(k) = \text{co} (Z(k))`.

* Update

    .. math::
        Z^{+}(k) = \left\{ z^{+}_l (k) \in Z(k) : z^{+}_l (k) \in \partial\mathcal{W}^{i+}(k) \right\};

  and find approximation subgradients :math:`\tilde{H}^{+}(k) = \{ \tilde{h}^{+}_1(k),...,\tilde{h}^{+}_{\tilde{L}^{+}}(k) \}` and constants :math:`C^{+}(k) = \{c^{+}_1(k),...,c^{+}_{\tilde{L}^+}(k) \}` such that

    .. math::
        co(Z^{+}(k)) = \bigcap_{l=1}^{\tilde{L}^{+}} \left\{ z: \tilde{h}_l^{+}(k) z \leq
        c^{+}_l (k) \right\},

  and :math:`\mathcal{W}^{i+} = \cup_{k} co(Z^{+}(k)) = \mathbf{B}^i (\mathcal{W}^i)`.

        


.. note ::
    
    As in the *outer approximation* methods, since :math:`A \times B` is a
    finite set of action profiles, we can evaluate the program :eq:`InnerBLP`
    as a special class of a nonlinear optimization problem--a *nonseparable*
    bilinear program (BLP)--for each fixed :math:`(a,b) \in A \times B`. Then
    we can maximize over the set :math:`A \times B`, by table look-up. Thus,
    the only difference computationally in the *inner approximation* method is
    the extra step of summarizing each inner step-correspondence :math:`\mathcal{W}^{i+}` by
    updates on:
    
        * approximation subgradients in each :math:`\tilde{H}^{+}(k)`;
    
        * levels in each :math:`C^{+}(k)`; and
    
        * vertices, :math:`Z^{+}(k)`,
    
    for every :math:`k \in \mathbf{K}`.

    
|



.. rubric:: Footnotes

.. [#footnote_subgradients_inner]

    Note that we have to let the *approximation subgradients* :math:`\tilde{H}` to possibly vary with domain partition elements :math:`Q_k`, as opposed to fixed *search subgradients* in :math:`H` used in the optimization step. This is because the former is endogenously determined by the extra convex hull operation taken to construct an inner step-correpondence :math:`\mathcal{W}^{i}` at each successive evaluation of the operator :math:`\mathbf{B}^i`.

.. topic:: Relevant functions :math:`\blacktriangleright`

    .. function:: Admit_Inner_LPset(self)

        Returns:

        * *Znew* :
    
            * A :math:`(L \times K)` numeric array containing elements
              :math:`z^{+}_l(k)` where :math:`k \in \mathbf{K}` and :math:`l \mapsto h_l \in
              H` are, respectively, a partition element of the correspondence
              domain :math:`D`, and, approximation subgradient in direction indexed
              by :math:`l \in \{1,...,L\}`.

    .. seealso:: **Punish_Inner**


|
|
|


.. [JYC2003] 
   Judd,  Kenneth L., Sevin Yeltekin and James Conklin (2003): "Computing Supergame Equilibria". *Econometrica*, 71(4), 1239-1254.

.. [SY2000] 
    Sleet, Christopher and Sevin Yeltekin (2000): "On the Computation of Value
    Correspondences". Unpublished. KGMS-MEDS, Northwestern University.

