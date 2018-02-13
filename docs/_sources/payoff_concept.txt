.. payoff_concept:

.. default-domain:: mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Approximating SSE operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The goal ahead is to *approximate* and provide computable representations of:

* each candidate correspondence :math:`\mathcal{W}: D \rightrightarrows
  \mathbb{R}^{\mathcal{Z}}`; and

* the operator :math:`\mathcal{W} \mapsto \mathbf{B}(\mathcal{W})`.

  
Conceptual
----------------


Recall :math:`\{Q_k\,|\,k=1,\ldots ,K\}` denotes a partition of :math:`D`, so
:math:`D=\bigcup_{k=1}^K Q_k`. An upper hemicontinuous, compact- and
convex-valued correspondence :math:`\mathcal{W}:D \rightrightarrows
\mathbb{R}^\mathcal{Z}` can be approximated by step-valued correspondences using the following procedures: Letting

.. math::
    :nowrap:

       \begin{equation*} 
            \omega^o_k (\lambda) := 
            \begin{cases}
            \text{co}\bigcup_{\lambda \in Q_k}\mathcal{W}(\lambda') & \text{if} \ (\lambda) \in Q_k,\\
            \emptyset & \text{otherwise},
            \end{cases}
        \end{equation*}

the correspondence defined by :math:`\mathcal{W}^o(\lambda) :=
\bigcup_{k}\omega^o_k (\lambda)` gives an *outer step-valued approximation*
of :math:`\mathcal{W}`.

Similarly, letting

.. math::
    :nowrap:

        \begin{equation*} 
            \omega^i_k (\lambda) := 
            \begin{cases}
             \bigcap_{\lambda \in Q_k}\mathcal{W}(\lambda) & \text{if} \ \lambda \in Q_k,
             \\
            \mathbb{R}^\mathcal{Z} & \text{otherwise},
            \end{cases}
        \end{equation*}

the correspondence defined by :math:`\mathcal{W}^i(\lambda) :=
\bigcap_{k}\omega^i_k (\lambda)` yields an *inner step-valued approximation*
of :math:`\mathcal{W}`. 

Practical
----------------


Since the convex-valued approximations
:math:`\mathcal{W}^o` and :math:`\mathcal{W}^i` are constant on each partition
element :math:`Q_k`, and there are :math:`K < +\infty` partition elements, these
approximations can be further approximated by constructing outer and inner
approximations for the sets :math:`\omega^o_k (\lambda)` and
:math:`\omega^i_k (\lambda)` using **convex polytopes**. 
Let :math:`\mathbb{S}^{N_z-1} := \left\{x \in \mathbb{R}^{N_z} : \| x \| = 1 \right\}` be the unit :math:`(N_z - 1)`-sphere where the norm :math:`\| \cdot \|` is given by :math:`\| x \|_{2} = \left(\sum_{j=1}^{N_z} x_{j}^2\right)^{1/2}`. Suppose we have finite sets of directional vectors: :math:`H := \{ h_l \in
\mathbb{S}^{N_z-1} : l = 1,...,L \}` and :math:`\tilde{H} := \{ \tilde{h}_l \in
\mathbb{S}^{N_z-1} : l = 1,...,L' \}`.
Let
:math:`\bar{\omega}^o_k (\lambda)` and :math:`\bar{\omega}^i_k (\lambda)`
denote the corresponding polytope approximations, respectively, of :math:`\omega^o_k (\lambda)` and
:math:`\omega^i_k (\lambda)`, where

.. math::
    :nowrap:

         \begin{equation*} 
            \bar{\omega}^o_k (\lambda) := 
            \begin{cases}
             \bigcap_{l=1}^{L}\{ z | h_l \cdot z \leq c_{l}^{o}(k) \} & \text{if} \ \lambda \in Q_k,\\
            \emptyset & \text{otherwise}
            \end{cases},
        \end{equation*}

and,

.. math::
    :nowrap:

         \begin{equation*} 
            \bar{\omega}^i_k (\lambda) := 
            \begin{cases}
             \bigcap_{l=1}^{L'}\{ z | \tilde{h}_l \cdot z \leq c_{l}^{i}(k) \} & \text{if} \ \lambda \in Q_k,\\
            \emptyset & \text{otherwise}
            \end{cases}.
        \end{equation*}


Let :math:`\bar{\mathcal{W}}^o := \cup_{k \in \mathbf{K}}\bar{\omega}^o_k`
and :math:`\bar{\mathcal{W}}^i :=  \cup_{k \in \mathbf{K}}\bar{\omega}^i_k` denote the resulting correspondences. One would like the "true" correspondence :math:`\mathcal{W}` to be "sandwiched" by polytope "step-correspondences" :math:`\bar{\mathcal{W}}^o` from the outside, and, by  :math:`\bar{\mathcal{W}}^i` from the inside. [#footnote_outerinner]_

.. math::
    :label: inner-outer-bound

            \bar{\mathcal{W}}^i \subset \mathcal{W}^i 
            \subset \mathcal{W} \subset \mathcal{W}^o
            \subset \bar{\mathcal{W}}^o.

The last statement :eq:`inner-outer-bound` is only true if the step-correspondence levels :math:`c_{l}^{o}(k)` and :math:`c_{l}^{i}(k)` are defined, respectively, as the maximal and minimal levels over each domain partition element :math:`Q_k`, in each direction :math:`h_l \in H` or :math:`\tilde{h}_l \in \tilde{H}`. [#footnote-ck]_

In the next section, we show how to construct these upper- and lower
bounding estimates :math:`c_{l}^{o}(k)` and :math:`c_{l}^{i}(k)` by using
stochastic global optimization programs and also *separable* bilinear program
formulations, when :math:`\mathcal{W}` represents a candidate guess of the
*symmetric sequential equilibrium* payoff correspondence in our class of
games.

.. rubric:: Footnotes

.. [#footnote-ck]
    In the context of our game, where :math:`\mathcal{W}` stands for a candidate guess of the equilibrium value correspondence, the last statement :eq:`inner-outer-bound` is only true if the step-correspondence levels :math:`c_{l}^{o}(k)` and :math:`c_{l}^{i}(k)` are defined, respectively, as the **globally** maximal and minimal values of each nonlinear programming problem (which is defined over each state-space partition element :math:`Q_k`, in each direction :math:`h_l \in H` or :math:`\tilde{h}_l \in \tilde{H}`) that summarizes the concept of *symmetric sequential equilibrium* of the game.


.. [#footnote_outerinner]
   This idea of providing both upper- and lower-bounding estimates of a given
   correspondence was first proposed by [JYC2003]_ in the computation of
   repeated games. Our proposed method is a modification of [SY2000]_
   who in turn extended [JYC2003]_ to the computation of value correspondences
   in dynamic games. Our
   contribution will be in the form of *bilinear programming formulations* as
   a practical and computable means of constructing these approximate correspondences.


