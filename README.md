QMC [![Build Status](https://img.shields.io/travis/matteosecli/QMC/BSc-Thesis.svg)](https://travis-ci.org/matteosecli/QMC) [![GitHub issues](https://img.shields.io/github/issues/matteosecli/QMC.svg)](https://github.com/matteosecli/QMC/issues) [![GitHub forks](https://img.shields.io/github/forks/matteosecli/QMC.svg)](https://github.com/matteosecli/QMC/network) [![GitHub stars](https://img.shields.io/github/stars/matteosecli/QMC.svg)](https://github.com/matteosecli/QMC/stargazers) [![GitHub license](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://raw.githubusercontent.com/matteosecli/QMC/BSc-Thesis/LICENSE) [![Unicorn Approval](http://img.shields.io/badge/unicorn-approved-ff69b4.svg)](https://www.youtube.com/watch?v=9auOCbH5Ns4) 
=======

About
-----
This is the repository for the QMC Project @ UiO, Autumn 2014. The current branch is an extension of the original project, meant to be presented as a BSc Thesis at the University of Trento

Introduction
------------
The aim of this project is to use the variational Monte Carlo (**VMC)** method to evaluate the ground state energy, one-body densities, expectation values of the kinetic and potential energies and single-particle energies of quantum dots with $N = 2$ and $N = 6$ electrons, so-called *closed shell systems*.
We will begin with two electrons confined in a pure two-dimensional isotropic harmonic oscillator potential, with and without the repulsive interaction, developing a simple trial wave-function and a basic program that runs our calculations. Then we are going to improve the trial wave-function adding the so-called *Jastrow factor*, and finally we are going to extend the entire class of wave-functions and the program itself to a given number of electrons (obviously, for a closed-shell system). After that we are going to improve our sampling, switching from a brute-force approach to the clever *importance sampling*.
In addition, we will present a method that implements analytical derivatives for calculating both the energy and the diffusive force of the importance sampling.	

Contacts
--------
Matteo Seclì (<secli.matteo@gmail.com>)