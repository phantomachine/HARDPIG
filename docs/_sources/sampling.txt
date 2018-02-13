.. sampling:

.. default-domain:: mat

Sampling on State-space Partitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

We will be drawing uniform sample vectors of the variable :math:`\lambda \in D` from each partition element :math:`Q_{k},k \in \mathbf{K}` of the domain :math:`D`. And then we will also use this in conjunction with the :ref:`xPoly`.

Each partition element in general is a tetrahedron simplex. We utilize a modified Markov-chain Monte Carlo sampler called the **Hit-and-Run Algorithm**, originally due to [Smi1984]_. This algorithm has a desirable property that it can (globally) reach any point in any arbitrarily given bounded set in :math:`\mathbb{R}^n` in one step. That is, there is a positive probability of sampling from any neighborhood in that set. Moreover, it is proven by [Lov1999]_ that the Hit-and-Run sampler converges fast (in probability) to a uniform distribution on convex bodies :math:`Q_k \subset \mathbb{R}^n`. [#footnote_hitandrun]_ [LV2003]_ note that this algorithm is the fastest in practice.


.. topic:: Hit-and-Run Algorithm

    Let :math:`S \subset \mathbb{R}^n` be a convex region that restricts sample
    realizations.
    The aim is to generate sample :math:`X := \{ s_m \}_{m=1}^{N_{sim}}` as a Markov Chain that is (asymptotically) uniformly distributed on :math:`S`. Define :math:`f(s)` by *any* continuous and strictly positive probability density function (pdf) on :math:`S`.

    * Start at a given point :math:`s` in the given set :math:`S`. Let
      :math:`m=1`.

    * Propose a new location :math:`s' = s + l d` by stepping away from :math:`s`
      according to a random direction-stepsize pair, :math:`(d,l)`, where the
      direction :math:`d` is uniformly distributed on the unit hypersphere
      :math:`\mathbb{S}^{n-1}` embedded in
      :math:`\mathbb{R}^n`; and the stepsize :math:`l \in \mathbb{R}` is drawn
      from a *proposal density* :math:`g_m(l | d, s)`.

    * Accept proposal move to :math:`s'` with *acceptance probability*
      :math:`\tilde{\alpha}_f(s,s')`, or, reject and stay at :math:`s` (i.e. set :math:`s' = s`) with probability
      :math:`1-\tilde{\alpha}_f(s,s')`.

    * Then set :math:`s_{m+1}` as :math:`s'`, and, repeat the procedure again from
      the first step, and let :math:`m = m + 1`.



To implement this simple algorithm, we need to define the functions :math:`\{ g_m(\cdot | d,s) : m \geq 1 \}`, and :math:`f` (which implies :math:`\tilde{\alpha}_f`) to ensure the necessary and sufficient (Kolmogorov) detailed balance condition holds (for the chain to be a *reversible Markov chain*):

        .. math ::
            g_m \left( \| s - s' \| ; \frac{s' - s}{\| s - s' \|}, s
            \right) \tilde{\alpha}_f(s,s')f(s) = 
            g_m \left( \| s - s' \| ; \frac{s - s'}{\| s - s' \|}, s'
            \right) \tilde{\alpha}_f(s',s)f(s').

This demands that the products of probabilities around every closed loop are the same in both directions around the loop.
 
* We can define :math:`f(s)` by *any* continuous and strictly positive
  probability density function (pdf) on :math:`S`.

* Let :math:`L_m := \left\{ l \in \mathbb{R} : s_m + l d_m \in S
  \right\}`. Define a conditional *proposal density* for each step :math:`m = 1,...,N_{sim}` by :math:`g_m ( l | d, s)`. 
  
    * Proposal densities that satisfy the *detailed balance condition* include the class of symmetric proposal density---i.e. :math:`g_m ( l | d, s) = \tilde{g}_m (l) = \tilde{g}(-l)` for all :math:`l \in \mathbb{R}`,  in which case
      
        .. math ::
            \tilde{\alpha}(s,s') = \min \left\{ \frac{f(s')}{f(s)} ,
            1\right\}.

    * Since we have also :math:`S` bounded, we can define a valid proposal
      density as
  
        .. math ::
            \tilde{g}_m (l) =  \frac{ \mathbf{1}_{ \{ l \in L_m \} } }{ \int_{\mathbb{R}} \mathbf{1}_{\{ u \in L_m \} } du }.

.. topic:: Example

    In our application, we will define a :math:`S := Q_k` for every :math:`k \in \mathbf{K}`.


See Section 6.3.1 of [KTB2011]_ for a generalized version of this simple algorithm.


.. topic:: Looking ahead

    Here we give a preview of the usage of uniform sampling from the convex partition elements. In :doc:`payoff_approximation` later, we show that in our class of dynamic games,
    the description of the *symmetric sequential equilibrium* operator (which is correspondence valued) involves solving many *non-separable* bilinear programs (BLP) of the form:

        .. math ::
            :nowrap:

            \begin{eqnarray*}
            \max_{\lambda,w} & u^{T} + \lambda^{T} Q w
            \\
            s.t. & \lambda \in \mathcal{F}_1(w),
            \\
                 & w \in \mathcal{F}_2(\lambda);
            \end{eqnarray*}

    where :math:`u \in  \mathbb{R}^{N_{z}}` is a vector of constants; :math:`\lambda,w \in \mathbb{R}^{N_{z}}` are the variables of interest; :math:`Q` is some :math:`(N_z \times N_z)` real matrix; and the constraint sets :math:`\mathcal{F}_1(w)` and :math:`\mathcal{F}_2(\lambda)` are convex polytopes which, respectively, depend on the choices of :math:`w` and :math:`\lambda`. [#footnote_nonseparable_bilinear]_ [#footnote_nonseparable_bilinear_source]_  

    We propose a Monte Carlo or stochastic approach to obtain :math:`\epsilon`-global (i.e. approximately global) optimization solutions to these non-separable bilinear programs. For now, notice that for each given realization of the random vector :math:`\lambda`, the nonseparable  BLP above can be reduced to standard linear programs (LP) in the variable :math:`w`. [#footnote_feasible_lambdas]_
   
.. rubric:: Footnotes

.. [#footnote_hitandrun] 
    [Lov1999]_ proves that the upper bound on the convergence rate is in
    polynomial time of :math:`\mathcal{O}(n^3)`.

.. [#footnote_nonseparable_bilinear] 
    A special (and textbook case) is where :math:`\mathcal{F}_1(w) \equiv \mathcal{F}_1` and :math:`\mathcal{F}_2(\lambda) \equiv \mathcal{F}_2`---i.e. each constraint set :math:`\mathcal{F}_1` and :math:`\mathcal{F}_2` do not vary, respectively, with the choice variables :math:`w` and :math:`\lambda`. This special case is known as a *separable bilinear* program, and, it nests quadratic programming as another special case. These problems are known to have a global solution--see [BM1993]_. Furthermore, successive approximation using branching-and-bounding methods--i.e. branching into subsets of the optimizer domain, then bounding the value function below by the solutions of linear programs on each subset of the function domain, and, above by the value from a local nonlinear optimizer--can be used to find the :math:`\epsilon`-global optimum: [McC1976]_ , [BM1993]_ and [HT1996]_

.. [#footnote_nonseparable_bilinear_source] 
    In the paper, we noted that in this class of games, the *source of
    bilinear nonseparability in the constraint sets* of :math:`\lambda` and :math:`w` is
    the utilitarian government's set of *incentive* or *promise-keeping* constraints.


.. [#footnote_feasible_lambdas]
    Note that by fixing each :math:`\lambda`, the constraint set :math:`\mathcal{F}_1(w)` will be redundant in the LP formulation within the stochastic global optimization scheme. Additionally, we will also require each realization :math:`\lambda` to be feasible according to some feasibility (e.g. a budget-balance) requirement(s): :math:`\lambda \in \mathcal{F}(Q_k) := \{ \lambda \in Q_k : \lambda b^{T} \geq 0, \forall b \in B \}`, where :math:`B` is some finite set of action profiles of the large (government) player.



.. topic:: Example (Sampling from :math:`Q_k` and :math:`Poly_{k'}`)

    The following figure (:ref:`figure_intersectsim`) shows an example of our usage of the Hit-and-Run algorithm in conjunction with our polytope intersection problems described earlier. For example, consider the (4,1)-panel in this figure. It shows the realizations of the random vectors :math:`\lambda P(a)`, where :math:`a = 4` denotes the fourth action profile in :math:`A`, that would end up in the various partition elements of :math:`(Q_3, Q_4, Q_5, Q_8)`, and given that each vector :math:`\lambda` is randomly drawn from the set :math:`Q_9`.


.. _figure_intersectsim:

.. figure:: _figures/intersectsim_a4.png
    :width: 100%
    :align: center
    :alt: Uniform samples of partitions
    :figclass: align-center

    Uniform samples and various polytope intersections :math:`(P(a_{4}), N_z = 3, K = 16)`




|
|
|
|
|

.. [BM1993] 
    Bennett, Kristin and Olvi L. Mangasarian (1993): "Bilinear Separation of Two Sets in n-Space". *Computational Optimization and Applications*, 2.

.. [HT1996]
    Horst, Reiner and Hoang Tuy (1996): "Global Optimization: Deterministic
    Approaches". Springer Verlag.

.. [KTB2011]
    Kroese, Dirk P., Thomas Taimre and Zdravko I. Botev (2011): *Handbook of
    Monte Carlo Methods*. Wiley Series in Probability and Statistics. Wiley.

.. [Lov1999] 
    Lovasz, Laszlo (1999): "Hit-and-run Mixes Fast". *Mathematical Programming*, Ser. A 86, 443-461.

.. [LV2003] 
    Lovasz, Laszlo and Santosh Vempala (2003): "Hit-and-Run is Fast and Fun",
    Technical Report, Microsoft Research, MSR-TR-2003-05.
    
.. [McC1976] 
    McCormick, Garth P. (1976): "Computability of Global Solutions to Factorable Nonconvex Programs: Part I. Convex underestimating problems". *Mathematical Programming*, 10, 147–175.

.. [Smi1984] 
    Smith, Robert L. (1986): "Efficient Monte-Carlo Procedures for Generating Points Uniformly Distributed over Bounded Regions". *Operations Research*, 32, 1296-1308. 
