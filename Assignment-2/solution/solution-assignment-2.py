from sys import exit, argv
import numpy as np


class component:

    def __init__(self, name):
        self.type = None
        self.name = name
        self.num_of_nodes = None
        self.n1 = None
        self.n2 = None
        self.voltage_source = None
        self.value = None
        self.voltage_source_n1 = None
        self.voltage_source_n2 = None

    def reversed_printing(self):
        if self.num_of_nodes == 2:
            print(comp.value + " " + comp.n2 + " " + comp.n1 + " " + comp.name)
        if self.num_of_nodes == 3:
            print(comp.value + " "+comp.voltage_source + " " +
                  comp.n2 + " " + comp.n1 + " " + comp.name)
        if self.num_of_nodes == 4:
            print(comp.value + " "+comp.voltage_source_n2 + " "+comp.voltage_source_n1 + " " +
                  comp.n2 + " " + comp.n1 + " " + comp.name)


def tokenizer(line: str):

    l_line = line.split(" ")
    names = {"R": "Resistor", "L": "Inductor", "C": "Capacitor", "V": "Independent Voltage source", "I": "Independent Current Source",
             "E": "Voltage Controlled Voltage Source", "G": "Voltage Controlled Voltage Source", "H": "Current Controlled Voltage Source", "F": "Current Controlled Current Source"}
    comp = component(l_line[0])
    comp.type = names[comp.name[0]]
    comp.num_of_nodes = len(l_line)-2
    if comp.num_of_nodes == 2:
        comp.n1 = l_line[1]
        comp.n2 = l_line[2]
        comp.value = l_line[3]
    if comp.num_of_nodes == 3:
        comp.n1 = l_line[1]
        comp.n2 = l_line[2]
        comp.voltage_source = l_line[3]
        comp.value = l_line[4]
    if comp.num_of_nodes == 4:
        comp.n1 = l_line[1]
        comp.n2 = l_line[2]
        comp.voltage_source_n1 = l_line[3]
        comp.voltage_source_n2 = l_line[4]
        comp.value = l_line[5]
    return comp


try:
    f = open(
        argv[1]
    )
    data = f.readlines()

    CIRCUIT = ".circuit"
    END = ".end"

    # to find the part where the first circuit definition starts and ends
    start = -1
    end = -2
    for i in range(len(data)):
        if CIRCUIT in data[i]:
            start = i+1
        if END in data[i]:
            end = i
            break

    if start >= end or start < 0:  # validating the circuit definition
        print("Invalid circuit definition")
        exit(0)

    lines = data[start:end]
    f.close()

    components = []
    for line in lines:  # tokenizing each component
        comp = tokenizer(line.split("#")[0].rstrip(" ").rstrip("\n"))
        components += [comp]

    # # printing each component in reversed order and in reversed format
    # for i in range(len(components)-1, -1, -1):
    #     comp = components[i]
    #     # print(comp.type)
    #     comp.reversed_printing()

    # solve the circuit


# following code runs if the filepath does not exist
except IOError:
    print("Invalid filename")
    exit()
