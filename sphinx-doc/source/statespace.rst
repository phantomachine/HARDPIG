.. statespace:

.. default-domain:: mat

State Space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Consider the *game state space* as :math:`D := \Delta (\mathcal{Z})`, the set of
all probability measures on the finite *individual state-space*
:math:`\mathcal{Z} := \{-N,...,-1,+1,...,+M\}`, where :math:`0<N,M <+\infty`.
Let :math:`N_{z}:=M+N \equiv |\mathcal{Z}|`.

===========================
Properties 
===========================

The set :math:`D` is:

#. is a unit simplex embedded in :math:`\mathbb{R}^{N_{z}}`:
    
.. math::
    :nowrap:
    
    \begin{equation*}
         \Delta (\mathcal{Z}) := \left\{ \lambda \in \mathbb{R}^{N_{z}}:
            \lambda_{i} \in [0,1], \forall i = 1,...,N_{z},
            \text{ and }
            \sum_{i=1}^{N_{z}} \lambda_{i} = 1
         \right\}
    \end{equation*}
    
#. represented by a convex polytope (i.e. a unit :math:`N_{z}`-simplex);

#. partitioned into :math:`K < +\infty` equal-area :math:`(N_{z}-1)`-simplices, :math:`Q_{k}, k \in \{ 1,...,K \} =: \mathbf{K}`.


.. topic:: Relevant functions :math:`\blacktriangleright`

    .. function:: simplex_tripart(K)

        Returns :math:`K` number of equal volume simplex partition elements of unit simplex :math:`D`, given by :math:`Q_{k}, k = 1,...,K`. 

The next figure--:ref:`figure_partition`--shows an example where :math:`K = 16` and :math:`N_z = 3`.

.. _figure_partition:

.. figure:: _figures/partition.png
    :width: 90%
    :align: center
    :alt: Partition Elements of D
    :figclass: align-center

    Example state-space partition scheme :math:`(N_z = 3, K = 16)`

   
