import numpy as np
import dimod
from dimod.reference.samplers import ExactSolver
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

# We create Q for the QUBO model and we also create J and h for the BinaryQuadraticModel
# h contains the diagonal values and J the triangular matix with the rest
Q = {}
J = {}
h = {}
values=[2,4,3,1]
values=[2,7,16,25,14,3,8,4,21,43,16,1]
c=0
for value in values:
  c+=value
i=0
j=0
while i < len(values):
    while j<len(values):
        if i==j:
            Q[(i,j)]=values[i]*(values[i]-c)
            h[i]=values[i]*(values[i]-c)
        else:
            Q[(i,j)]=values[i]*values[j]
            J[(i,j)]=values[i]*values[j]
        j+=1
    i+=1
    j=0
i=0
j=0


#We can create the model from the qubo matrix
model=dimod.BinaryQuadraticModel.from_qubo(Q)
print("The model we are going to solve in QUBO is")
print(model)
print()
#Or directly using h and J
model = dimod.BinaryQuadraticModel(h, J, 0.0, dimod.BINARY)
print("The model we are going to solve in BinaryQuadratic is")
print(model)
print()


# First we solve using the exact solver

sampler = ExactSolver()
solution = sampler.sample(model)
print("The exact solution is ")
print(solution)
print()


# Now with *simulated annealing*

sampler = dimod.SimulatedAnnealingSampler()
response = sampler.sample(model, num_reads=10)
print("La solucion con simulated annealing es")
print(response)
print()


# Finally we solve it again using the *quantum annealer* selecting explicitly the computer to use


sampler = EmbeddingComposite(DWaveSampler(solver='Advantage_system1.1'))
sampler_name = sampler.properties['child_properties']['chip_id']
response = sampler.sample(model, num_reads=5000)
print("The solution with D-Wave's quantum annealer called ",sampler_name,"en")
print(response)
print()

# The same with the other *annealer*

sampler = EmbeddingComposite(DWaveSampler(solver='DW_2000Q_6'))
sampler_name = sampler.properties['child_properties']['chip_id']
response = sampler.sample(model, num_reads=5000)
print("The solution with D-Wave's quantum annealer called ",sampler_name,"en")
print(response)
print()
