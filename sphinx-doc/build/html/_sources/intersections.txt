.. intersections:

.. default-domain:: mat


==============================================
Intersections with State-space Partitions
==============================================

For every :math:`k \in \mathbf{K}` and its associated simplicial partition element :math:`Q_{k} \subset D` with
positive volume, the
set-valued image :math:`P(a)(Q_{k})`:

#. is another :math:`N_{z}`-simplex contained in the unit :math:`N_{z}`-simplex :math:`D`; and

#. intersects with:

    * at least one partition element :math:`Q_{k'}` where :math:`k' \in \mathbf{K}` and

    * at most all partition elements :math:`Q_{1}, ..., Q_{K}`;

.. _xPoly:

---------------------------------
Polytope intersection problems
---------------------------------

Denote

.. math ::
         \mathbf{I}(a,k) := \left\{ k' \in \mathbf{K} :
                                    P(a)(Q_{k}) \cap Q_{k'} \neq \emptyset
                            \right\}, 
                            \qquad
                            \forall a \in A, k \in \mathbf{K},    

as the sets of indexes to respective partition-elements---i.e. :math:`k'
\mapsto Q_{k'}`---that contain non-empty
intersections with each simplicial image :math:`P(a)(Q_{k})`. Each nonempty
intersection, induced by each :math:`(a,k) \in A \times \mathbf{K}` and :math:`P(a)`, is described by
    
.. math ::
    Poly_{k'(a,k)} := \left\{ \lambda' \in D :
                                 \lambda' \in P(a)(Q_{k}) \cap Q_{k'}
                                 \neq \emptyset
                                 \text{ and }
                                 k' \in \mathbf{I}(a,k) 
                         \right\}.  

.. note ::
    Each intersection :math:`Poly_{k'}`, for each
    :math:`k \in \mathbf{K}` and each :math:`a \in A`, is a
    *polytope*, and is *at least a simplex*, and is a subset of partition
    element :math:`Q_{k'}`, where :math:`k' := k'(a,k)`.

These nonempty intersections are such that

.. math ::
    \bigcup_{k' \in \mathbf{I}(a,k)} Poly_{k'} :=: P(a)(Q_{k}).

.. topic:: Example
    
    If :math:`N_{z} = 3`, then :math:`D` is a unit 2-simplex, and
    each :math:`Poly_{k'}` can be a polygon or a triangular subset in :math:`D`.




