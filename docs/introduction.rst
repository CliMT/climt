.. highlight:: python

=============
Introduction
=============

climt (pronounced *klimt*) is an attempt to build a climate modelling
infrastructure that is easy to use, easy to understand and easy to learn.

Most climate model components are written in fortran for performance reasons.
For that very reason, it is difficult to change model configurations and 
behaviour easily, which is something scientists tend to do all the time during
their research. The earth-science community is converging towards Python as the
language of choice for data analysis tasks, thanks to Python's flexibility and
emphasis on clean, readable code. climt aims to use Python for climate modelling
for these very reasons -- clearly documented components, self documenting
model scripts, and a flexible configuration system will make climate modelling more
reproducible and make the learning curve less steep.

climt is aimed at a wide spectrum of users -- from students who are curious to learn
about the climate system to researchers who want state-of-the-art components. climt
aims to provide multiple levels of abstraction which will allow the user to tradeoff
ease of use vs. flexibility, depending on their particular needs and experience in 
modelling and Python programming.
