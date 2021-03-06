import matplotlib.pyplot as plt
import numpy as np


data = np.loadtxt("test 1.tsv")
calc_sep_vecs = np.transpose(data[:,0])
calc_forces = np.transpose(data[:,1])
write_state = np.transpose(data[:,2])
calc_dist = np.transpose(data[:,3])
update_state = np.transpose(data[:,4])

data = np.loadtxt("test 2.tsv")
calc_sep_vecs_fast = np.transpose(data[:,0])
calc_forces_fast = np.transpose(data[:,1])
write_state_fast = np.transpose(data[:,2])
calc_dist_fast = np.transpose(data[:,3])

data = np.loadtxt("test 3.tsv")
calc_sep_vecs_faster = np.transpose(data[:,0])
calc_forces_faster = np.transpose(data[:,1])
write_state_faster = np.transpose(data[:,2])
calc_dist_faster = np.transpose(data[:,3])

plt.loglog(calc_sep_vecs)
plt.loglog(calc_sep_vecs_fast)
plt.loglog(calc_sep_vecs_faster)
plt.title("calc_sep_vecs() function")
plt.xlabel("Particle Number")
plt.ylabel("Time Taken (s)")
plt.grid(which='both')
plt.tight_layout()
plt.savefig("Figures/calc_sep_vecs.pdf")
plt.show()

plt.loglog(calc_forces)
plt.loglog(calc_forces_fast)
plt.loglog(calc_forces_faster)
plt.title("calc_forces() function")
plt.xlabel("Particle Number")
plt.ylabel("Time Taken (s)")
plt.grid(which='both')
plt.tight_layout()
plt.savefig("Figures/calc_forces.pdf")
plt.show()

plt.loglog(write_state)
plt.loglog(write_state_fast)
plt.loglog(write_state_faster)
plt.title("write_state() function")
plt.xlabel("Particle Number")
plt.ylabel("Time Taken (s)")
plt.grid(which='both')
plt.tight_layout()
plt.savefig("Figures/write_state.pdf")
plt.show()

plt.loglog(calc_dist)
plt.loglog(calc_dist_fast)
plt.loglog(calc_dist_faster)
plt.title("calc_dist() function")
plt.xlabel("Particle Number")
plt.ylabel("Time Taken (s)")
plt.grid(which='both')
plt.tight_layout()
plt.savefig("Figures/calc_dist.pdf")
plt.show()
'''
plt.plot(update_state)
plt.title("update_state() function")
plt.xlabel("Particle Number")
plt.ylabel("Time Taken (s)")
plt.tight_layout()
plt.savefig("Figures/update_state.pdf")
plt.show()
'''
'''
data2 = np.loadtxt("energies.tsv")
kin = np.transpose(data2[:,0])
pot = np.transpose(data2[:,1])
tot = np.transpose(data2[:,2])

plt.plot(kin,label="Kinetic")
plt.plot(pot,label="Potential")
plt.plot(tot,label="Total")
plt.xlabel("Timestep")
plt.ylabel("Energy")
plt.tight_layout()
plt.legend()
plt.savefig("Figures/Energy.pdf")
plt.show()
'''