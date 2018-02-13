.. tute1:

%%%%%%%%%%%%%%%%%%%%
Tutorial 1
%%%%%%%%%%%%%%%%%%%%

Begin by browsing through appendices A and B of `Econometric Theory <http://johnstachurski.net/personal/emet_book.html>`_.
(Feel free to skim with the
intention of going back to it again.)
Step through the code you find there and type it into the R interpreter.
Try variations of what you see:  

* If you see ``x <- 1``, try ``x <- 2``
* If you see ``col="red"``, try ``col="blue"``

If variations lead to errors,

* Read the error message and try to understand what the problem is
* Try `googling <http://www.google.com>`_ the error message (cut and paste into google)

Once you've finished (or feel like it), try the following exercises.  If you
get stuck, either go back and read those appendices again, or, if you must, have a
look at the solutions below.


Exercises
===========

Exercise 1
------------

In **one** line of code, create a list called ``rnums`` with two elements, where

* the names of the elements are ``a`` and ``b`` respectively,
* the first element is a vector of 1,000 draws from a standard normal distribution,
* the second element is a vector of 1,000 draws from a uniform distribution on [0,1]

You should now find that ``mean(rnums$b)`` returns a value close to 0.5.  Why?

Exercise 2
------------

Next, in **one** line of code, create a histogram from all 2,000 random variables in ``rnums``.


Exercise 3
------------

A Gaussian random walk is a sequence of random variables given by

.. math::
    X_t = \sum_{i=1}^t W_i  \quad t = 1, 2, 3, \ldots
    :label: grwalk



where each ``W_i`` is an (independent) draw from the standard normal distribution.

* Investigate the function ``cumsum()`` by typing ``?cumsum`` at the prompt
* Using this function, in one line of code, generate a vector ``x`` of length 1,000 containing a simulated Gaussian random walk
* Plot the random walk as a time series using a blue line
 
Exercise 4
------------

Here's a high-school maths problem:  

* Fred is born on the first of January, 1970
* Sally is born on the first of January, 1992

In which year will Fred be exactly 3 times as old as Sally?
How old will they be?

A bit of simple maths determines that 

* Year = 2003
* Fred's age = 33
* Sally's age = 11

Your task is to write a computer program that outputs the same results.
Write it in a script file.

If you are new to programming, the problem is not trivial.  Try. 
If you can't do it then wait for the solution to be posted, and study it.

Here's what I want you to do:

First, write a function ``freds_age()`` such that, given integer ``year`` greater
than 1970, ``freds_age(year)`` returns the age of Fred in that year.

Second, write a function ``sallys_age()`` such that, given integer ``year`` greater
than 1992, ``sallys_age(year)`` returns the age of Sally in that year.

Test these functions, and make sure they work as expected.

Now, create a variable ``year`` with the value of 1993 (when Sally is 1).

Next, write a ``while`` loop that tests whether Fred's age is > 3 times Sally's age
in the current ``year``.  If yes, ``year`` is incremented by one.  If not, the loop terminates.

Finally, write three lines to print out the value of ``year`` and Fred and Sally's ages
in that year.



Solutions
===========

Solution to Exercise 1
------------------------



.. code-block:: r

    rnums <- list(a=rnorm(1000), b=runif(1000))




Solution to Exercise 2
------------------------



.. code-block:: r

    hist(c(rnums$a, rnums$b))




Solution to Exercise 3
------------------------



.. code-block:: r

    x <- cumsum(rnorm(1000))
    plot(x, col="blue", type="l")

Solution to Exercise 4
------------------------




.. code-block:: r

    sally_born <- 1992
    fred_born <- 1970
    
    sallys_age <- function(y) return(y - sally_born)
    freds_age <- function(y) return(y - fred_born)
    
    year <- 1993
    while (3 * sallys_age(year) < freds_age(year)) {
        year <- year + 1
    }
    
    cat("Year: ", year, "\n") 
    cat("Sally's age: ", sallys_age(year), "\n") 
    cat("Fred's age: ", freds_age(year), "\n") 


Ted's MATLAB solution
-------------------------

.. code-block:: matlab

	if strcmp(plot_option, 'none') ~= 1          
            simplex2dset_partdraw(D, T_new, dt, ic, DT, K_DT, K, J, n);
            fprintf('... DONE plotting, K = %i elements\n', 4^n)
        end
