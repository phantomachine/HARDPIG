.. actions_transitions:

.. default-domain:: mat

Actions and Transitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Let :math:`\tilde{A} \ni a_{j}` denote the finite action set of each individual (small)
player :math:`j \in \mathcal{Z}`. Then :math:`A := \tilde{A}^{\mathcal{Z}} \ni a`, denotes the
:math:`N_{z}`-copies of an individual's action space---i.e. the set of action
profiles :math:`a`. The finite  set :math:`A` has :math:`N_{a} := |A|` number of
action profiles.

For each action profile :math:`a \in A`, its associated Markov transition
probability function is a linear operator :math:`P(a): D \rightarrow D`.


.. topic:: Relevant functions :math:`\blacktriangleright`

    .. function:: A_ProfileSet(self)
 
        Returns the finite set :math:`A` of :math:`N_{a}` action profiles
        of individual small players. A numeric array of size :math:`N_{a} \times N_{z}`.
        
    .. function:: TransProbA(self)

        Returns the finite set of Markov maps :math:`P(a): D \rightarrow D`, one
        for every :math:`a \in A`. A numeric array of size :math:`N_{z} \times
        N_{z} \times N_{a}`.




