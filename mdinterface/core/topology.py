#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 14:05:39 2024

@author: roncofaber
"""

import copy
#%%

class Topology(object):
    
    def __init__(self, resname=None, formula=None):
        
        self._id = None
        self._formula = formula
        self._resname = resname
        
        return
    
    @property
    def id(self):
        return self._id
    
    @property
    def formula(self):
        return self._formula
    
    @property
    def resname(self):
        return self._resname
    
    def set_id(self, value):
        self._id = value
        return
    
    def set_formula(self, value):
        self._formula = value
        
    def set_resname(self, value):
        self._resname = value
        
    def copy(self):
        return copy.deepcopy(self)

#%%

class Atom(Topology):
    def __init__(self, symbol, label=None, eps=None, sig=None):
        
        super(Atom, self).__init__()
        
        self.symbol  = symbol
        
        self.eps = eps
        self.sig = sig
        
        if label is None:
            label = symbol

        self.set_label(label)
        
        return
    
    def set_label(self, value):
        self._label = value
        return
    
    @property
    def label(self):
        return self._label
    
    @property
    def extended_label(self):
        return self._label + "_" + self.resname
    
    def __repr__(self):
        return f"{self.__class__.__name__}({self.label}, eps={self.eps}, sig={self.sig})"
    
    # TODO CHANGE BACK IF DOESN'T WORK
    def __eq__(self, other):
        return self.eps == other.eps and self.sig == other.sig and self.symbol == other.symbol

    # def __eq__(self, other):
    #     return self.eps == other.eps and self.sig == other.sig and self.label == other.label
    
#%%

class Bond(Topology):
    def __init__(self, a1, a2, kr=None, r0=None):
        super(Bond, self).__init__()
        self._a1 = a1
        self._a2 = a2
        self.kr = kr
        self.r0 = r0
        return

    @property
    def symbols(self):
        return self._a1, self._a2

    def __repr__(self):
        return f"{self.__class__.__name__}({'-'.join(self.symbols)},kr={self.kr}, r0={self.r0})"
    
    @property
    def elements(self):
        return self._a1.split('_')[0], self._a2.split('_')[0]

    def __eq__(self, other):
        return (self.kr == other.kr and self.r0 == other.r0 and
                (self.elements == other.elements or self.elements == other.elements[::-1]))


class Angle(Topology):
    def __init__(self, a1, a2, a3, kr=None, theta0=None):
        super(Angle, self).__init__()
        self._a1 = a1
        self._a2 = a2
        self._a3 = a3
        self.kr = kr
        self.theta0 = theta0
        return

    @property
    def symbols(self):
        return self._a1, self._a2, self._a3
    
    @property
    def elements(self):
        return self._a1.split('_')[0], self._a2.split('_')[0], self._a3.split('_')[0]
    
    def __repr__(self):
        return f"{self.__class__.__name__}({'-'.join(self.symbols)}, kr={self.kr}, theta0={self.theta0})"

    def __eq__(self, other):
        return (self.kr == other.kr and self.theta0 == other.theta0 and
                (self.elements == other.elements or self.elements == other.elements[::-1]))


class Dihedral(Topology):
    def __init__(self, a1, a2, a3, a4, A1=None, A2=None, A3=None, A4=None, A5=None):
        super(Dihedral, self).__init__()
        self._a1 = a1
        self._a2 = a2
        self._a3 = a3
        self._a4 = a4
        self._values = [A1, A2, A3, A4, A5]
        return

    @property
    def symbols(self):
        return self._a1, self._a2, self._a3, self._a4

    @property
    def values(self):
        return self._values
    
    @property
    def elements(self):
        return self._a1.split('_')[0], self._a2.split('_')[0],\
            self._a3.split('_')[0], self._a4.split('_')[0]

    def __repr__(self):
        symbols = '-'.join(self.symbols)
        values = '|'.join([str(ii) for ii in self._values])
        return f"{self.__class__.__name__}({symbols}, A={values})"

    def __eq__(self, other):
        return (self.values == other.values and
                (self.elements == other.elements or self.elements == other.elements[::-1]))

    def write(self, fout):
        atype = "{}-{}-{}-{}".format(*self.symbols)
        if self._values[-1] is not None:
            value = "{:>7.4f}  {:>7.4f}  {:>7.4f}  {:>7.4f}  {:>7.4f}".format(*self.values)
        else:
            value = "{:>7.4f}  {:>7.4f}  {:>7.4f}  {:>7.4f}".format(*self.values[:-1])
        fout.write("{:>5}    {}  #  {:<8} | {}\n".format(self.id, value, atype, self.resname))


class Improper(Topology):  # cvff improper style
    def __init__(self, a1, a2="x", a3="x", a4="x", K=None, d=None, n=None):
        super(Improper, self).__init__()
        
        self._a1 = a1
        self._a2 = a2
        self._a3 = a3
        self._a4 = a4
        
        self._K = K
        
        assert d in [1, -1], "Wrong 'd' value"
        self._d = d
        
        assert n in [0,1,2,3,4,5,6], "Wrong 'n' value"
        self._n = n
        
        return

    @property
    def symbols(self):
        return self._a1, self._a2, self._a3, self._a4

    @property
    def values(self):
        return self._K, self._d, self._n
    
    def __repr__(self):
        symbols = '-'.join(self.symbols)
        values = '|'.join([str(ii) for ii in self.values])
        return f"{self.__class__.__name__}({symbols}, K,d,n={values})"

    def __eq__(self, other):
        return (self.values == other.values and
                self._a1.split('_')[0] == other._a1.split('_')[0])




