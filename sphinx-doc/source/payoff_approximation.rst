.. payoff_approximation:

.. default-domain:: mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Equilibrium Payoff Correspondence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

In :doc:`../statespace` we constructed a partition of the simplex :math:`D:= \Delta(\mathcal{Z})`. Now, we let :math:`D` be the domain of the equilibrium payoff correspondence.
The task ahead is to approximate the equilibrium value correspondence :math:`\mathcal{V}: D
\rightrightarrows \mathbb{R}^{\mathcal{Z}}` using
convex-valued step correspondences. 

----------------------------------------------
Background
----------------------------------------------

Notation reminder:

* Action profile of small players on :math:`[0,1]`, :math:`a \in A`. (Assume :math:`A` is a
  finite set.) Each small player takes on a *personal state* :math:`j \in
  \mathcal{Z} := \{ -N,...,-1,+1,...,+M\}` at each date :math:`t \in
  \mathbb{N}`.

* Actions of large player (:math:`G`), :math:`b \in B`. :math:`B := \{b \in
  \mathbb{R}^{\mathcal{Z}}: -m \leq b(j) \leq \bar{m}, \text{ for } j > 0,
  \text{ and, } 0 \leq b(j) \leq \bar{m}, \text{ for } j < 0, \forall j \in
  \mathcal{Z}\}` is a
  finite set and contains vectors :math:`b` that are physically feasible
  (but not necessarily government-budget feasible in all states).

* Extended payoff vector space, :math:`\mathbb{R}^{\overline{\mathcal{Z}}}`,
  where :math:`\overline{\mathcal{Z}} := \mathcal{Z} \cup \{G\}`.
  
* Probability distribution of small players on finite set
  :math:`\mathcal{Z}`, :math:`\lambda \in D := \Delta(\mathcal{Z})`.

* Profile of continuation values of agents, :math:`w \in
  \mathbb{R}^{\mathcal{Z}}`.

* Transition probability matrix at action profile :math:`a`, :math:`P(a)`
    
* Individual :math:`j \in \mathcal{Z}`, given action :math:`a(j)` faces transition probability distribution, :math:`p^j(a(j)) \in P(a)`

* Flow payoff profile, :math:`v_j(a,b):= u(c^b(j))-\phi(a_j)`, where
  
    * :math:`v(a,b):=(v_j(a,b))_{j\in \mathcal{Z}}`;

    * Utility-of-consumption function, :math:`u(\cdot)`; and

    * Disutility-of-effort/action function, :math:`-\phi(\cdot)`.

* Public date-:math:`t` history, :math:`h^t := \{ \lambda^t, x^t, b^{t-1} \}_{t \geq 0}`, where

    * :math:`\lambda^t = (\lambda_0,...,\lambda_t)` is a history of agent
      distributions up to and include that of date :math:`t`;

    * :math:`x^t = (x_0,...,x_t)` where :math:`x_t` is a date-:math:`t` realization of the
      random variable :math:`X_t \sim_{i.i.d.} \mathbf{U}([0,1])`; and
      
    * :math:`b^{t-1} = (b_0,...,b_{t-1})` is a history of government policy
      actions up to the end of date :math:`t-1` and let :math:`\{b_{-1}\} =
      \emptyset`.

  Also, :math:`h^0 := (\lambda_0, x_0)`.

.. _def-consistency:

.. topic:: Definition (Consistency)   

    Let :math:`\mathcal{W}:D \rightrightarrows \mathbb{R}^{\overline{\mathcal{Z}}}` be a compact- and convex-valued correspondence having the property that :math:`w(G)= \sum_{j\in\mathcal{Z}}\lambda(j)w(j)` for all :math:`(\lambda,w)\in\text{graph}(\mathcal{W})`. A vector :math:`(b,a,\lambda',w)\in B \times A \times\Delta(\mathcal{Z}) \times\mathbb{R}^{\overline{\mathcal{Z}}}` is **consistent with respect to** :math:`\mathcal{W}` at :math:`\lambda` if
    
    #. :math:`-\lambda b^{T} \geq 0`;

    #. :math:`\lambda'=\lambda P(a)`;

    #. :math:`w\in\mathcal{W}(\lambda')`; and
    
    #. For all :math:`j\in\mathcal{Z}`, :math:`a(j)\in\text{argmax}_{a'}\left\{(1-\delta) \left[u(c^b(j))-\phi(a')\right]+\delta\mathbb{E}_{p^j(a')}[w(i)]\right \}`.

.. _def_admissibility:

.. topic:: Definition (Admissibility)

    For :math:`(\lambda,b) \in D \times B` let 
        .. math::
            \pi(\lambda,b):=\min_{(a',\lambda'',w')} 
            \left[(1-\delta)\sum_{j\in\mathcal{Z}}\lambda(j)[u(c^b(j))-\phi(a'(j))]
            + \delta\sum_{j\in\mathcal{Z}}\lambda''(j)w'(j)\right],
            
    subject to :math:`(b,a',\lambda'',w')` is consistent with respect to :math:`\mathcal{W}(\lambda)`. Let :math:`(\tilde{a}(\lambda,b),\tilde{\lambda}'(\lambda,b),\tilde{w}(\lambda,b))` denote the solutions to the corresponding minimization problem. A vector :math:`(b,a,\lambda',w)\in B \times A \times D \times\mathbb{R}^{\overline{\mathcal{Z}}}` is said to be **admissible with respect to** :math:`\mathcal{W}(\lambda)` if

    #. :math:`(b,a,\lambda',w)` is consistent with respect to :math:`\mathcal{W}(\lambda)`; and
    
    #. :math:`(1-\delta)\sum_{j\in\mathcal{Z}}\lambda(j)[u(c^b(j))-\phi(a(j))]\delta\sum_{j\in\mathcal{Z}}\lambda'(j)w(j)\geq \max_{b'\in B(\lambda)}\pi(\lambda,b')`, where :math:`B(\lambda) := \{ b \in B :- \lambda b^{T} \geq 0 \}`.


.. topic:: Admissible payoff vectors

    The payoff vector defined by an admissible vector :math:`(b,a,\lambda',w)` at :math:`\lambda` is given by

    .. math::
        E_G(b,a,\lambda',w)(\lambda)&=(1-\delta)\sum_{j\in\mathcal{Z}}\lambda(j)[u(c^b(j))-\phi(a(j))]+ \delta\sum_{j\in\mathcal{Z}}\lambda'(j)w(j)
        \\
        E_j(b,a,\lambda',w)(\lambda)&=(1-\delta) \left[u(c^b(j))-\phi(a(j))\right] + \delta \mathbb{E}_{p^j(a(j))}[w(i)].

    Note that :math:`E_G(b,a,\lambda',w)(\lambda)=\sum_{j\in\mathcal{Z}}\lambda(j)E_j(b,a,\lambda',w)(\lambda)`. 
    
In the paper, we proved the following:

.. topic:: SSE Recursive Operator

    A SSE is a strategy profile :math:`\sigma \equiv \{ \alpha_t, \beta_t\}_{t \geq 0}` such that given initial game state :math:`\lambda_0`, for all dates :math:`t \geq 0`, and all public histories :math:`h^t := \{ \lambda^t, x^t, b^{t-1} \}_{t \geq 0}`, :math:`a := \alpha_t(h^t, b_t)` and :math:`b := \beta_t(h^t)`, and, if :math:`\mathcal{V}` is the SSE payoff correspondence, then :math:`\mathcal{V}` is the largest fixed point that satisfies the recursive operator
    
    .. math::
        \mathbf{B}(\mathcal{V})(\lambda):=\text{co}\{E(b,a,\lambda',w)(\lambda)\,\vert\,
        (b,a,\lambda',w) &\text{ is admissible w.r.t.} \mathcal{V}(\lambda)\},    
    
    where :math:`\text{co}` denotes the convex hull of a set.

    The object of interest can be found recursively: :math:`\mathcal{V} = \lim_{n
    \rightarrow +\infty} \mathbf{B}^n(\mathcal{W}_0)`, for any initial
    convex-valued and compact correspondence :math:`\mathcal{W}_0`.

.. note ::
    Given state :math:`\lambda` and agent payoff vector
    :math:`w \in \mathbb{R}^\mathcal{Z}` determine a unique corresponding
    government payoff given by :math:`\lambda \cdot w`. We can thus ignore the
    government payoff when defining the equilibrium value correspondences and
    their approximations, and restrict their codomain to
    :math:`\mathbb{R}^\mathcal{Z}`.


