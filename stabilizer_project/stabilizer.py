"Contains the classes and function to manipulate stabilizer and graph states"
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.providers.aer import QasmSimulator
from qiskit.visualization import plot_histogram
from qiskit.quantum_info import StabilizerState, Pauli
import matplotlib
import matplotlib.pyplot as plt

class Stab:
    '''
    This is a class that encodes the stabilizer state in terms of its stabilizers
    
    :param size: Number of qubits
    :type size: int

    :param stabs: The stabilizers, either in a string or a list, in the format 'XX,-YY' or '[XX,-YY]' (case sensitive). Optional, defaults to 'XX,ZZ'
    :type stabs: list or string, optional
    '''
    def __init__(self, n = 2, stabs = 'XX,ZZ'):
        """Constructor method

        """
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
        """
        Returns the size of the stabilizer (the number of qubits)

        :return: The size of the stabilizer
        :rtype: int
        """
        return self.size
    def tableau(self):
        """
        Converts the stabilizers to a tableau and signvector

        :return: A list contained the tableau and the signvector
        :rtype: list
        """
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
        """
        Resets the stabilizer and new tableau associated with it

        :param newstabs: The new stabilizers
        :type newstabs: string
        """
        try:
            self.stab = newstabs.split(',')
        except:
            self.stab = newstabs
        list = self.tableau()
        self.tab = list[0]
        self.signvector = list[1]
    def clifford(self,type,q1,q2=None):
        """
        Applies a clifford gate to the stabilizer

        :param type: The clifford gate to be operated, 'H', 'X', 'Y', 'Z', 'CNOT', 'CZ', or 'S'
        :type type: string

        :param q1: The qubit to operate on, or the control qubit for entangling gates
        :type q1: int

        :param q2: The qubit to target, defaults to None
        :type q2: int
        """
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
        """
        Prints the tableau and the signvector

        
        """
        print(self.tab)
        print(self.signvector)
    def gaussian(self):
        """
        Reorder the stabilizers to make the circuit builder easier

        
        """
        pass
    def circuit_builder(self):
        """
        Uses reverse operations to build the stabilizer state

        :return: A list of operations to take a standard state to the given stabilizer state
        :rtype: list        
        """
        size = self.size
        rev_operations = []
        reference = np.copy(self.tab)
        for i in range(self.size):
            if self.tab[i,i]==1:
                continue
            elif self.tab[i,self.size+i]==1:
                self.clifford('h',i)
                rev_operations.append('h,'+str(i))
        for i in range(self.size):
            for j in range(self.size):
                if i != j:
                    if self.tab[i,j]==1:
                        self.clifford('CNOT',i,j)
                        rev_operations.append('CNOT,'+str(i)+','+str(j))
        for i in range(self.size):
            if self.tab[i,self.size+i]==1:
                if self.signvector[i]==0:
                    rev_operations.append('s,'+str(i))
                    self.clifford('s',i)
                    self.clifford('z',i)
                elif self.signvector[i]==1:
                    rev_operations.append('sdag,'+str(i))
                    self.clifford('s',i)
        print(rev_operations)
        print(self.tab)
        self.tab = np.copy(reference)
        operations = rev_operations.reverse()
        return operations


def grapher(edgelist):
    """
    Function that can graph a graph state provided an edgelist

    Parameters
    ----------
    edgelist : list
        A list that denotes all the connections in a graph state. The edgelist is a nested list, with each inner list containing two elements, the numbers of the qubits that are connected.

    Returns
    -------
    circuit : QuantumCircuit
        A Qiskit quantum circuit that encodes the circuit that generates that graph
    """
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

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    grapher([[0,1],[1,2],[2,3],[3,4],[4,5],[5,0]])
