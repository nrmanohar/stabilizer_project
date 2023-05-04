"Contains the classes and function to manipulate stabilizer and graph states"
import numpy as np
from qiskit import QuantumCircuit, ClassicalRegister
from qiskit.quantum_info import StabilizerState

class Stabilizer:
    '''
    This is a class that encodes the stabilizer state in terms of its stabilizers. If no input is given, it will initialize a bell state. If only the n is given, it will initialize n qubits in the 0 state
    
    :param n: Number of qubits
    :type n: int, Optional

    :param stabs: The stabilizers, either in a string or a list, in the format 'XX,-YY' or '[XX,-YY]' (case sensitive). Optional
    :type stabs: list or string, optional

    :ivar size: The number of qubits, initial value: n
    :ivar __stabs: The stabilizers of the state, initial value: stabs (note, this is a dunder attribute, can't be directly called outside the class. There's a method to do that instead)
    :ivar tab: The tableau of the state
    :ivar signvector: The signvector of the state
    :ivar gauss: A nxn Gaussian matrix (used for empty_column calculations)
    
    '''
    def __init__(self, n = None, stabs = None, edgelist = None):
        """Constructor method

        """
        if edgelist is None:    
            if n is None and stabs is None:
                n = 2
                stabs = 'XX,ZZ'
            
            elif n is not None and stabs is None:
                stabs = []
                for i in range(n):
                    str = ''
                    for j in range(n):
                        if i==j:
                            str = str+'Z'
                        else:
                            str = str+'I'
                    stabs.append(str)
            

            self.size = n
            try:
                self.__stab = stabs.split(',')
            except:
                self.__stab = stabs
            list = self.tableau()
            self.tab = list[0]
            self.signvector = list[1]
            while not self.square():
                print("Invalid input, number of qubits not equal to number of stabilizers")
                n = int(input("Number of Qubits "))
                stabs = input("Stabilizers ")
                self.size = n
                try:
                    self.__stab = stabs.split(',')
                except:
                    self.__stab = stabs
                list = self.tableau()
                self.tab = list[0]
                self.signvector = list[1]
            while self.empty_column():
                print("Invalid input, free qubit (all stabilizers for some qubit is the identity)")
                n = int(input("Number of Qubits "))
                stabs = input("Stabilizers ")
                self.size = n
                try:
                    self.__stab = stabs.split(',')
                except:
                    self.__stab = stabs
                list = self.tableau()
                self.tab = list[0]
                self.signvector = list[1]
            while not self.commuter():
                print("Invalid Inputs, Stabilizers do not commute")
                n = int(input("Number of Qubits "))
                stabs = input("Stabilizers ")
                self.size = n
                try:
                    self.__stab = stabs.split(',')
                except:
                    self.__stab = stabs
                list = self.tableau()
                self.tab = list[0]
                self.signvector = list[1]
            while not self.linear_independence():
                print("Invalid Inputs, Stabilizers are not independant")
                n = int(input("Number of Qubits "))
                stabs = input("Stabilizers ")
                self.size = n
                try:
                    self.__stab = stabs.split(',')
                except:
                    self.__stab = stabs
                list = self.tableau()
                self.tab = list[0]
                self.signvector = list[1]
        else:
            self.graph_state()

    def square(self):
        toggler = True
        for i in range(len(self.__stab)):
            str = self.__stab[i]
            str = str.lstrip('-')
            if len(str)!=len(self.__stab):
                return False
        return toggler

    def commuter(self):
        """
        Tests whether the stabilizers commute with each other

        :return: Whether or not they commute
        :rtype: boolean
        """
        for i in range(self.size):
            toggler=0
            for j in range(i+1,self.size):
                for k in range(self.size):
                    if self.tab[i,k]==self.tab[j,k] and self.tab[i,k+self.size]==self.tab[j,k+self.size]:
                        toggler = toggler
                    elif (self.tab[i,k+self.size]==0 and self.tab[i,k]==0) or (self.tab[j,k+self.size]==0 and self.tab[j,k]==0):
                        toggler = toggler
                    else:
                        toggler = toggler+1
                if toggler%2 != 0:
                    return False
        return True
    def num_qubits(self):
        """
        Returns the size of the stabilizer (the number of qubits)

        :return: The size of the stabilizer
        :rtype: int
        """
        return self.size
    
    def graph_state(self,edgelist = [[0,1],[1,2],[2,3],[3,4],[4,0]]):
        num=0
        for i in range(len(edgelist)):
            for j in range(len(edgelist[i])):
                if edgelist[i][j]>num:
                    num = edgelist[i][j]
        self.size = num+1
        tab = np.zeros(2*self.size*self.size)
        tab = tab.reshape(self.size,2*self.size)
        for i in range(self.size):
            tab[i,i]=1
        for i in range(len(edgelist)):
            q1 = edgelist[i][0]
            q2 = edgelist[i][1]
            tab[q1,q2+self.size]=1
            tab[q2,q1+self.size]=1
        sign = np.zeros(self.size)
        self.tab = tab
        self.signvector = sign
    
    def tableau(self):
        """
        Converts the stabilizers to a tableau and signvector

        :return: A list contained the tableau and the signvector
        :rtype: list
        """
        tab = np.zeros(2*self.size*self.size)
        tab = tab.reshape(self.size,2*self.size)
        sign = np.zeros(self.size)
        for i in range(len(self.__stab)):
            if self.__stab[i][0]=='-':
                sign[i]=1
                self.__stab[i]=self.__stab[i][1:]
            for j in range(len(self.__stab[i])):
                if self.__stab[i][j]=='I':
                    pass
                elif self.__stab[i][j]=='X':
                    tab[i,j]=1
                elif self.__stab[i][j]=='Z':
                    tab[i,j+self.size]=1
                elif self.__stab[i][j]=='Y':
                    tab[i,j]=1
                    tab[i,j+self.size]=1
                else:
                    print('Invalid Stabilizer')
        return [tab,sign]
    def stabilizers(self):
        """
        Returns a list of the stabilizers of the state, as per the tableau

        :return: A list of operations to take a standard state to the given stabilizer state
        :rtype: list  
        
        """
        self.__stab = []
        for i in range(self.size):
            str = ""
            if self.signvector[i]==1:
                str = str+"-"
            for j in range(self.size):
                if self.tab[i,j]==0 and self.tab[i,j+self.size]==0:
                    str = str+"I"
                elif self.tab[i,j]==1 and self.tab[i,j+self.size]==0:
                    str = str+"X"
                if self.tab[i,j]==0 and self.tab[i,j+self.size]==1:
                    str = str+"Z"
                if self.tab[i,j]==1 and self.tab[i,j+self.size]==1:
                    str = str+"Y"
            self.__stab.append(str)
        return self.__stab
    def new_stab(self,size=None,newstabs=None):
        """
        Resets the stabilizer and new tableau associated with it

        :param size: The size of the new state
        :type size: int (optional)

        :param newstabs: The new stabilizers
        :type newstabs: string or list
        """
        if size is None and newstabs is None:
            size = 2
            newstabs = 'XX,ZZ'
        
        if size is not None and newstabs is None:
            newstabs = []
            for i in range(size):
                str = ''
                for j in range(size):
                    if i==j:
                        str = str+'Z'
                    else:
                        str = str+'I'
                newstabs.append(str)
        self.size = size
        try:
            self.__stab = newstabs.split(',')
        except:
            self.__stab = newstabs
        list = self.tableau()
        self.tab = list[0]
        self.signvector = list[1]
        while not self.square():
            print("Invalid input, number of qubits not equal to number of stabilizers")
            n = int(input("Number of Qubits "))
            stabs = input("Stabilizers ")
            self.size = n
            try:
                self.__stab = stabs.split(',')
            except:
                self.__stab = stabs
            list = self.tableau()
            self.tab = list[0]
            self.signvector = list[1]
        while self.empty_column():
            print("Invalid input, free qubit (all stabilizers for some qubit is the identity)")
            n = int(input("Number of Qubits "))
            stabs = input("Stabilizers ")
            self.size = n
            try:
                self.__stab = stabs.split(',')
            except:
                self.__stab = stabs
            list = self.tableau()
            self.tab = list[0]
            self.signvector = list[1]
        while not self.commuter():
            print("Invalid Inputs, Stabilizers do not commute")
            n = int(input("Number of Qubits "))
            stabs = input("Stabilizers ")
            self.size = n
            try:
                self.__stab = stabs.split(',')
            except:
                self.__stab = stabs
            list = self.tableau()
            self.tab = list[0]
            self.signvector = list[1]
        while not self.linear_independence():
            print("Invalid Inputs, Stabilizers are not independant")
            n = int(input("Number of Qubits "))
            stabs = input("Stabilizers ")
            self.size = n
            try:
                self.__stab = stabs.split(',')
            except:
                self.__stab = stabs
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
                self.clifford('h',q2)
                self.clifford('cnot',q1,q2)
                self.clifford('h', q2)
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
        Generates an array that contains information about where stabilizers are known
        
        """
        self.gauss = np.zeros(self.size*self.size)
        self.gauss = self.gauss.reshape(self.size,self.size)
        for i in range(self.size):
            for j in range(self.size):
                if self.tab[i,j]==1 or self.tab[i,j+self.size]==1:
                    self.gauss[i,j]=1
    
    def empty_column(self):
        """
        Tests whether there are any empty stabilizers (free qubits)

        :return: Whether there is an empty column or not
        :rtype: boolean
        """
        self.gaussian()
        zed = self.gauss.sum(axis=0)
        empty = False
        for i in range(self.size):
            if zed[i]==0:
                empty = True
        return empty

    def linear_independence(self):
        """
        Checks if the generators are linearly independent

        """
        rank = np.linalg.matrix_rank(self.tab)
        rank = int(rank)
        if rank == self.size:
            return True
        else:
            return False

    def row_add(self,row1,row2):
        """
        Multiplies two stabilizers in the tableau together, specifying a new stabilizer, and puts them into the second row

        """
        stabs=self.stabilizers()
        stab1 = stabs[row1]
        stab2 = stabs[row2]
        stab1 = stab1.lstrip('-')
        stab2 = stab2.lstrip('-')
        phase = np.ones(self.size, dtype=complex)
        for i in range(self.size):
            if stab1[i]=='Z' and stab2[i]=="X":
                phase[i]=np.complex(1j)*phase[i]
            elif stab1[i]=='X' and stab2[i]=="Z":
                phase[i]=np.complex(-1j)*phase[i]
            elif stab1[i]=='Y' and stab2[i]=="Z":
                phase[i]=np.complex(1j)*phase[i]
            elif stab1[i]=='Z' and stab2[i]=="Y":
                phase[i]=np.complex(-1j)*phase[i]
            elif stab1[i]=='X' and stab2[i]=="Y":
                phase[i]=np.complex(1j)*phase[i]
            elif stab1[i]=='Y' and stab2[i]=="X":
                phase[i]=np.complex(-1j)*phase[i]
        toggler = 1
        for i in range(self.size):
            toggler = toggler*phase[i]
        toggler = (1-1*np.real(toggler))/2
        self.signvector[row2]=(self.signvector[row1]+self.signvector[row2]+toggler)%2
        for i in range(2*self.size):
            self.tab[row2,i]=(self.tab[row2,i]+self.tab[row1,i])%2
     
            

    def circuit_builder(self):
        """
        Uses reverse operations to build the stabilizer state

        :return: A Qiskit circuit that makes the stabilizer
        :rtype: QuantumCircuit        
        """
        reference = np.copy(self.tab)
        sign = np.copy(self.signvector)
        rev_operations = []

        broken = False

        for i in range(self.size):
            if self.tab[i,i]==0:
                if self.tab[i,i+self.size]==1:
                    rev_operations.append(['H',i])
                    self.clifford('H',i)
            if self.tab[i,i]==0:
                for j in range(i+1,self.size):
                    if self.tab[j,i]==1:
                        self.tab[[i,j]]=self.tab[[j,i]]
                        self.signvector[[i,j]]=self.signvector[[j,i]]
                        break
            if self.tab[i,i]==0:
                for j in range(i+1,self.size):
                    if self.tab[j,i+self.size]==1:
                        self.tab[[i,j]]=self.tab[[j,i]]
                        self.signvector[[i,j]]=self.signvector[[j,i]]
                        rev_operations.append(['H',i])
                        self.clifford('H',i)
                        break
            if self.tab[i,i]==0:
                for j in range(i):
                    if self.tab[j,i+self.size]==1:
                        self.row_add(j,i)
                        rev_operations.append(['H',i])
                        self.clifford('H',i)
                        break
            if self.tab[i,i]==0:
                broken = True
                break
            elif self.tab[i,i]==1:
                for j in range(self.size):
                    if self.tab[i,j]==1 and j!=i:
                        rev_operations.append(["CNOT",i,j])
                        self.clifford("CNOT",i,j)
                

        if broken:
            self.tab = np.copy(reference)
            self.signvector = np.copy(sign)
            print("Something went wrong in the building procedure. Check your stabilizers and maybe reformat them and try again")
            return None

        for i in range(self.size):
            if self.tab[i,i+self.size]==1:
                rev_operations.append(["S",i])
                self.clifford("S",i)

        for i in range(self.size):
            for j in range(self.size):
                if self.tab[i,j+self.size]==1:
                    rev_operations.append(["CZ",i,j])
                    self.clifford("CZ",i,j)

        for i in range(self.size):
            self.clifford('H',i)
            rev_operations.append(['H',i])
        
        for i in range(self.size):
            if self.signvector[i]==1:
                rev_operations.append(['X',i])
                self.clifford('X',i)
        
        self.tab = np.copy(reference)
        self.signvector = np.copy(sign)
        
        rev_operations.reverse()

        circuit = QuantumCircuit(self.size)

        
        for i in range(len(rev_operations)):
            if rev_operations[i][0]=='H':
                circuit.h(rev_operations[i][1])
            elif rev_operations[i][0]=='S':
                circuit.s(rev_operations[i][1])
                circuit.z(rev_operations[i][1])
            elif rev_operations[i][0]=='X':
                circuit.x(rev_operations[i][1])
            elif rev_operations[i][0]=='Y':
                circuit.y(rev_operations[i][1])
            elif rev_operations[i][0]=='Z':
                circuit.z(rev_operations[i][1])
            elif rev_operations[i][0]=='CNOT':
                circuit.cnot(rev_operations[i][1],rev_operations[i][2])
            elif rev_operations[i][0]=='CZ':
                circuit.cz(rev_operations[i][1],rev_operations[i][2])
        return circuit

    def draw_circuit(self, style = 'mpl', save = None):
        """
        Draws a circuit that can generate the given stabilizer state (requires matplotlib and pylatexenc package)

        :param style: The type of output, 'mpl' for matplotlib, 'text' for ASCII drawing, 'latex_source' for raw latex output
        :type style: String, optional. Defaults to 'mpl'

        :param save: If you want to save the file to something (optional)
        :type save: String

        """
        if style == 'mpl':
            try:
                import matplotlib
                import matplotlib.pyplot as plt
                try:
                    circ = self.circuit_builder()
                    circ.draw(output = style, filename = save)
                    plt.show()
                except:
                    print("pylatexenc likely not installed")
            except:
                print("matplotlib package not installed")
        elif style == 'text':
            circ = self.circuit_builder()
            circ.draw(output = style, filename = save)
        elif style == 'latex_source':
            circ = self.circuit_builder()
            circ.draw(output = style, filename = save)

    
    def qiskit_stabilizers(self):
        """
        Asks Qiskit to return the stabilizers

        :return: A qiskit stabilizer state representation
        :rtype: StabilizerState (qiskit)

        """
        circ = self.circuit_builder()
        stab = StabilizerState(circ)
        return stab

    def stabilizer_measurement(self):
        """
        A circuit to implement the circuit and then to measure the associated stabilizers.

        :return: A qiskit circuit for measureing stabilizer
        :rtype: StabilizerState (qiskit)

        """
        qs = QuantumCircuit(2*self.size)
        bits = []
        for i in range(self.size):
            bits.append(i)
        reg = ClassicalRegister(self.size)
        qs.add_register(reg)
        stabs = self.stabilizers()
        for i in range(self.size):
            for j in range(self.size):
                if stabs[i][j]=='X':
                    qs.cx(self.size+i,j)
                elif stabs[i][j]=='Z':
                    qs.cz(self.size+i,j)
                elif stabs[i][j]=='Y':
                    qs.cz(self.size+i,j)
        for i in range(self.size):
            if self.signvector[i]==1:
                qs.x(self.size+i)
        for i in range(self.size):
            qs.measure(i+self.size,i)

        return qs
    
    def build_and_measure(self):
        circ = self.circuit_builder()
        qs = QuantumCircuit(2*self.size)
        bits = []
        for i in range(self.size):
            bits.append(i)
        qs = qs.compose(circ,bits)
        reg = ClassicalRegister(self.size)
        qs.add_register(reg)
        qs.barrier()
        stabs = self.stabilizers()
        for i in range(self.size):
            for j in range(self.size):
                if stabs[i][j]=='X':
                    qs.cx(self.size+i,j)
                elif stabs[i][j]=='Z':
                    qs.cz(self.size+i,j)
                elif stabs[i][j]=='Y':
                    qs.cz(self.size+i,j)
        for i in range(self.size):
            if self.signvector[i]==1:
                qs.x(self.size+i)
        for i in range(self.size):
            qs.measure(i+self.size,i)
        
        return qs
        
def grapher(edgelist):
    """
    Function that can graph a graph state provided an edgelist

    Parameters
    ----------
    edgelist : list
        A list that denotes all the connections inp a graph state. The edgelist is a nested list, with each inner list containing two elements, the numbers of the qubits that are connected.

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
