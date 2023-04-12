import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.providers.aer import QasmSimulator
from qiskit.visualization import plot_histogram
from qiskit.quantum_info import StabilizerState, Pauli
import matplotlib
import matplotlib.pyplot as plt

class Stab:
    '''
    This is a simple class to calculate the hamiltonian for a spin configuration
    :param size: Number of qubits
    :type size: int, optional
    :param stabs: The stabilizers, either in a string or a list
    :type mu: string or list, optional
    '''
    def __init__(self, n = 2, stabs = 'XX,ZZ'):
        self.size = n
        try:
            self.stab = stabs.split(',')
        except:
            self.stab = stabs
        self.signvector = np.zeros(n)
        list = self.tableau()
        self.tab = list[0]
        self.signvector = list[1]
    def num_qubits(self):
        return self.size
    def tableau(self):
        tab = np.zeros(2*self.size*self.size)
        tab = tab.reshape(self.size,2*self.size)
        sign = np.zeros(self.size)
        for i in range(len(self.stab)):
            if self.stab[i][0]=='-':
                sign[i]=1
                self.stab[i]=self.stab[i][1:]
            for j in range(len(self.stab[i])):
                if self.stab[i][j]=='I':
                    pass
                elif self.stab[i][j]=='X':
                    tab[i,j]=1
                elif self.stab[i][j]=='Z':
                    tab[i,j+self.size]=1
                elif self.stab[i][j]=='Y':
                    tab[i,j]=1
                    tab[i,j+self.size]=1
                else:
                    print('Invalid Stabilizer')
        return [tab,sign]
    def new_stab(self,newstabs):
        try:
            self.stab = newstabs.split(',')
        except:
            self.stab = newstabs
        list = self.tableau()
        self.tab = list[0]
        self.signvector = list[1]
    def clifford(self,type,q1,q2=None):
        if type.lower() == 'h':
            for i in range(self.size):
                alpha = self.tab[i,q1]
                beta = self.tab[i,q1+self.size]
                self.tab[i,q1]=beta
                self.tab[i,q1+self.size]=alpha
                self.signvector[i] = (self.signvector[i]+(alpha*beta))%2
        elif type.lower() == 'cnot':
            if q2 == None:
                print('Recall method and specify second qubit')
            else:
                for i in range(self.size):
                    self.tab[i,q2] = (self.tab[i,q1]+self.tab[i,q2])%2
                    self.tab[i,self.size+q1] = (self.tab[i,q1+self.size]+self.tab[i,q2+self.size])%2
                    if self.tab[i,q1]==1 and self.tab[i,q2+self.size]==1:
                        if self.tab[i,q2]==self.tab[i,self.size+q1]:
                            self.signvector[i]=(self.signvector[i]+1)%2
        elif type.lower() == 'z':
            for i in range(self.size):
                if self.tab[i,q1]==1:
                    self.signvector[i]=(self.signvector[i]+1)%2
        elif type.lower() == 'x':
            for i in range(self.size):
                if self.tab[i,q1+self.size]==1:
                    self.signvector[i]=(self.signvector[i]+1)%2
        elif type.lower() == 'y':
            for i in range(self.size):
                if (self.tab[i,q1]+self.tab[i,q1+self.size])==1:
                    self.signvector[i]=(self.signvector[i]+1)%2
        elif type.lower() == 's':
            for i in range(self.size):
                if self.tab[i,q1]==1:
                    self.signvector[i]=(self.signvector[i]+self.tab[i,q1+self.size])%2
                    self.tab[i,q1+self.size] = (self.tab[i,q1+self.size]+1)%2
        elif type.lower() == 'cz':
            if q2 == None:
                print('Recall method and specify second qubit')
            else:
                self.tab[:, [q2+self.size, q2]] = self.tab[:, [q2, q2+self.size]]
                for i in range(self.size):
                    if self.tab[i,q2]==1:
                        if self.tab[i,self.size+q2]==1:
                            self.signvector[i] = (self.signvector[i]+1)%2
                for i in range(self.size):
                    self.tab[i,q2] = (self.tab[i,q1]+self.tab[i,q2])%2
                    self.tab[i,self.size+q1] = (self.tab[i,q1+self.size]+self.tab[i,q2+self.size])%2
                    if self.tab[i,q1]==1 and self.tab[i,q2+self.size]==1:
                        if self.tab[i,q2]==self.tab[i,self.size+q1]:
                            self.signvector[i]=(self.signvector[i]+1)%2
                self.tab[:, [q2+self.size, q2]] = self.tab[:, [q2, q2+self.size]]
                for i in range(self.size):
                    if self.tab[i,q2]==1:
                        if self.tab[i,self.size+q2]==1:
                            self.signvector[i] = (self.signvector[i]+1)%2
        else:
            print("Something went wrong, make sure you inputted a valid type. Valid types are 'H' for Hadamard, 'S' for the phase gate, 'CNOT' for the Control Not, 'CZ' for the Control Z.")
    def report(self):
        print(self.tab)
        print(self.signvector)



def grapher(edgelist):
    num=0
    for i in range(len(edgelist)):
        for j in range(len(edgelist[i])):
            if edgelist[i][j]>num:
                num = edgelist[i][j]
    circuit = QuantumCircuit(num+1, num+1)
    for i in range(num+1):
        circuit.h(i)
    for i in range(len(edgelist)):
        circuit.cz(edgelist[i][0],edgelist[i][1])
    stab = StabilizerState(circuit)
    print(stab)
    return circuit


def canvas(with_attribution=True):
    """
    Placeholder function to show example docstring (NumPy format).

    Replace this function and doc string for your own project.

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from.

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution.
    """

    quote = "The code is but a canvas to our imagination."
    if with_attribution:
        quote += "\n\t- Adapted from Henry David Thoreau"
    return quote


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())
