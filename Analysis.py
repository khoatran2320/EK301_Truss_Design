import numpy as np
from scipy import linalg
from scipy.io import savemat, loadmat
import math


# Coded by Khoa Tran for

# Truss design preliminary analysis 

# Truss specifications:
#     joint to joint span:        6.5in to 15in
#     truss span:                 28in
#     load to pin support span:   13in
#     total virtual cost:         less than $275

##### Member properties #####
#
#   F(l) = 3127.25 / l^2        - where F(l) is the maximum compressive force prior to buckling (oz)
#                               - l is the length of the member (in)
#                               - the uncertainty is 1.01 oz
#                               
#   weight of material: .0843 oz/in


# Cost of truss:
#     Cost = c1j + c2l    where c1 = $10/joint, c2 = $1/in, j = total number of joints, and l = total length of all members

class Analysis:
    def __init__(self):
        self.cost_per_joint = 10
        self.mem_cost_per_in = 1
        self.num_mems = 0
        self.num_joints = 0
        self.total_mem_len = 0
        self.buckling_const = 3127.25
        self.weight_per_oz = .0843
        self.uncertainty = 1.01
        self.c_mat = None     #connection matrix
        self.sx_mat = None    #support matrix for x reaction forces
        self.sy_mat = None    #support matrix fro y reaction forces
        self.x_vect = None    #joint location x vector
        self.y_vect = None    #joint location y vector
        self.l_vect = None    #load vector
        self.a_mat = None     #coefficient matrix
        self.t_vect = None    #tension matrix solution
        self.truss_cost = None 
        self.r_vect = None
        self.load = None
        self.len_vect = None
        self.p_crit_vect = None
        self.weight_vect = None
        self.dead_load_vect = None
        self.failure_vect = None
        self.crit_member = None
        self.max_load = None
        self.load_cost_ratio = None
        self.account_for_dl = False
        self.dead_w_v = None
        self.cal_weak = False
        self.cal_strong = False

##### Calculations #####
    def calc_dead_load_vect(self):
        self.dead_load_vect = []
        for _ in range(2*self.num_joints):
            self.dead_load_vect.append(0)
        # print(self.weight_vect)
        for i in range(self.num_joints):
            joint_load = 0
            members_in_joint = self.find_mems_in_joint(i)
            # print(i, members_in_joint)
            for member in members_in_joint:
                if member == -1:
                    continue
                joint_load -= self.weight_vect[member]/2
            self.dead_load_vect[i+self.num_joints] = joint_load
        for i in range(2*self.num_joints):
            if self.l_vect[i] != 0:
                self.dead_load_vect[i] += float(self.l_vect[i])
        self.dead_w_v = [float(i) for i in self.dead_load_vect]
        try:
            l = np.array(self.dead_load_vect)
            l = l.reshape(2*self.num_joints,1)
            a = np.array(self.a_mat)
            self.dead_load_vect = np.linalg.solve(a, l)
        except:
            print("Unable to solve, exiting")
            exit(127)
        # print(self.dead_load_vect)
    def calc_failure_vect(self):
        self.calc_r_vect()
        self.failure_vect = []
        r_vect = self.r_vect
        t_vect = self.dead_load_vect if self.account_for_dl else self.t_vect

        for i in range(self.num_mems):
            if t_vect[i] < 1 and self.r_vect[i] != 0:
                if self.cal_weak:
                    self.failure_vect.append((self.p_crit_vect[i] - self.uncertainty)/r_vect[i])
                elif self.cal_strong:
                    self.failure_vect.append((self.p_crit_vect[i] + self.uncertainty)/r_vect[i])
                else:
                    self.failure_vect.append(self.p_crit_vect[i]/r_vect[i])
            else:
                self.failure_vect.append(0)
        c_vect = [abs(i) if i < 0 else 1.7976931348623157e+308 for i in self.failure_vect ]

        self.max_load = min(c_vect)
        self.crit_member = c_vect.index(min(c_vect))
        self.load_cost_ratio = self.max_load / self.truss_cost

    def calc_buckling_f(self, len):
        return self.buckling_const / (len**2)

    def calc_joint_dist(self, j1, j2):
        return math.sqrt( ((j2[1] - j1[1])**2) + ((j2[0] - j1[0])**2))     
    
    def find_mems_in_joint(self, joint):
        mem = []
        for i in range(self.num_mems):
            if self.c_mat[joint][i] == 1:
                mem.append(i)
            else:
                mem.append(-1)
        return mem
    
    def find_joint_pos_for_member(self, member):
        first_joint_found = False 
        j1 = []
        j2 = []
        for i in range(self.num_joints):
            if self.c_mat[i][member] == 1:
                if first_joint_found:
                    j2 = [self.x_vect[0][i], self.y_vect[0][i]]
                else:
                    first_joint_found = True
                    j1 = [self.x_vect[0][i], self.y_vect[0][i]]
        return j1, j2

    def find_joint_pos(self, joint):
        return [self.x_vect[0][joint], self.y_vect[0][joint]]

    def construct_A_mat(self):
        A = []
        for i in range(2*self.num_joints):
            A.append([0])

        for i in range(self.num_joints):
            arrx = []
            arry = []
            joint_pos = self.find_joint_pos(i)
            members_in_joint = self.find_mems_in_joint(i)
            for member in members_in_joint:
                if member == -1:
                    arrx.append(0)
                    arry.append(0)
                    continue
                j1_pos, j2_pos = self.find_joint_pos_for_member(member)
                if joint_pos == j2_pos:
                    coefx = (j2_pos[0] - j1_pos[0]) / self.calc_joint_dist(j1_pos, j2_pos)
                    coefy = (j2_pos[1] - j1_pos[1]) / self.calc_joint_dist(j1_pos, j2_pos)
                    arrx.append(coefx)
                    arry.append(coefy)
                elif joint_pos == j1_pos:
                    coefx = (j1_pos[0] - j2_pos[0]) / self.calc_joint_dist(j1_pos, j2_pos)
                    coefy = (j1_pos[1] - j2_pos[1]) / self.calc_joint_dist(j1_pos, j2_pos)
                    arrx.append(coefx)
                    arry.append(coefy)
                else:
                    arrx.append(0)
                    arry.append(0)
            for el in self.sx_mat[i]:
                arrx.append(el)
            for el in self.sy_mat[i]:
                arry.append(el)
            
            A[i] = arrx
            A[i+self.num_joints] = arry
        self.a_mat = A

    def solve(self):
        try:
            self.construct_A_mat()
            l = self.l_vect
            a = np.array(self.a_mat)
            self.t_vect = np.linalg.solve(a, l)
        except:
            print("Unable to solve, the matrix might be singular")
            
            # try:
            #     pinv = np.linalg.pinv(a)
            #     self.t_vect = pinv.dot(l)
            # except:
            #     print("Unable to solve, exiting")
            exit(127)

    def calc_truss_cost(self):
        self.truss_cost = self.cost_per_joint * self.num_joints + self.mem_cost_per_in * self.total_mem_len

    def calc_r_vect(self):
        if self.account_for_dl:
            self.r_vect = [float(t / self.load) for t in self.dead_load_vect]
        else:
            self.r_vect = [float(t / self.load) for t in self.t_vect]

##### Utilities #####
    def print_formatted_per_mem(self):
        # if not self.t_vect:
        self.solve()
        # if not self.truss_cost:
        self.calc_truss_cost()
        if self.account_for_dl:
            self.calc_dead_load_vect()
        self.calc_failure_vect()
        out = []
        t_v = self.dead_load_vect if self.account_for_dl else self.t_vect
        for el in t_v:
            out.append(el[0])
        out = [float(s) for s in out]
        out = [[s*-1, "C"] if s < 0 else [0,""] if s == 0 else [s, "T"] for s in out]
        self.construct_A_mat()
        l = [ -1*self.max_load if i != 0 else 0 for i in self.l_vect]
        
        if self.account_for_dl:
            for i in range(2*self.num_joints):
                if abs(self.dead_w_v[i]) > abs(self.load):
                    l[i] += -1*abs(abs(self.dead_w_v[i]) - abs(self.load))
                else:
                    l[i] += self.dead_w_v[i]
        # print(l)
        a = np.array(self.a_mat)
        t = np.linalg.solve(a, l)
        # print(t)
        print("Member number, Member length (in.), Tension (T) or Compression (C), Buckling strength (oz.) and uncertainty (%), Force at max load (oz.)")
        p_out = []
        for i in range(self.num_mems):
            p_out = []
            p_out.append(f"m{i+1}")
            p_out.append(f"{self.len_vect[i]:.2f}")
            p_out.append(out[i][1])
            p_out.append(f"{self.p_crit_vect[i]:.2f} with {self.uncertainty/self.p_crit_vect[i]*100:.2f}% uncertainty")
            p_out.append(f"{abs(t[i]):.2f}")
            for el in p_out[:-1]:
                print(el, end=", ")
            print(p_out[-1])
        
    def mat_save(self, filename = "TrussDesign1_KhoaDarinKevin_A1.mat"):
        to_save = {'C': self.c_mat, 'Sx': self.sx_mat, 'Sy': self.sy_mat, 'X': self.x_vect, 'Y': self.y_vect, 'L': self.l_vect}
        savemat(filename, to_save)
    
    def mat_load(self, filename = "TrussDesign1_KhoaDarinKevin_A1.mat"):
        res = loadmat(filename)
        try:
            self.c_mat = res['C']
            self.sx_mat = res['Sx']
            self.sy_mat = res['Sy']
            self.x_vect = res['X']
            self.y_vect = res['Y']
            self.l_vect = res['L']
        except:
            print("Unable to parse file")
    def print_output(self):
        if self.t_vect == None:
            self.solve()
        if self.truss_cost == None:
            self.calc_truss_cost()
        if self.account_for_dl:
            self.calc_dead_load_vect()
        self.calc_failure_vect()
        out = []
        t_vect = self.dead_load_vect if self.account_for_dl else self.t_vect
        for el in t_vect:
            out.append(el[0])
        out = [float(s) for s in out]
        out = [[s*-1, "C"] if s < 0 else [0,""] if s == 0 else [s, "T"] for s in out]
        
        
        print("\\% EK301, Section A1, Group1: Khoa T., Darin S., Kevin P., 4/1/2021.")
        print(f"Load: {self.load} oz")
        print("Member forces in oz")
        for i in range(self.num_mems):
            print(f"m{i+1}: %.3f ({out[i][1]})" %out[i][0] if out[i][0] > 1e-9 else f"m{i+1}: 0")
        print("Reaction forces in oz:")
        print(f"Sx1: %.3f" % out[self.num_mems][0])
        print(f"Sy1: %.3f" % out[self.num_mems+1][0])
        print(f"Sy2: %.3f" % out[self.num_mems+2][0])
        print("Cost of truss: $%.2f" %self.truss_cost)
        print("Theoretical max load/cost ratio in oz/$: %.5f" % self.load_cost_ratio)
        print(f"Critical member: m{self.crit_member+1} with critical buckling strength {self.p_crit_vect[self.crit_member]:.2f} oz and uncertainty {self.uncertainty/self.p_crit_vect[self.crit_member]*100:.2f}%")
        print(f"Max load: {self.max_load:.2f} oz with uncertainty {(self.uncertainty/self.p_crit_vect[self.crit_member])*self.max_load:.4f} oz")

    def print_truss_details(self):
        print("C matrix: ")
        print(self.c_mat)
        print("\nSx matrix: ")
        print(self.sx_mat)
        print("\nSy matrix: ")
        print(self.sy_mat)
        print("\nX vector: ")
        print(self.x_vect)
        print("\nY vector: ")
        print(self.y_vect)
        print("\nL vector: ")
        print(self.l_vect)
        print("\nA matrix:\n")
        for i in range(self.num_joints*2):
            for j in range(self.num_mems+3):
                print("%.3f\t" % self.a_mat[i][j], end='')
            print("\n")
        print("\n")

        print("T vector:\n")
        print(self.t_vect)
        print("\nNumber of members: ", self.num_mems)
        print("Number of joints: ", self.num_joints)
        

##### Checks #####
    def check_valid_truss(self):
        if not (self.num_mems == 2*self.num_joints - 3):
            print("Invalid truss!!")
            exit(127)
    
    def check_valid_c_mat(self):
        temp = self.c_mat.sum(axis=0)
        if not np.equal(2, temp)[0]:
            print("Invalid connection matrix!!")
            exit(127)

    def check_valid_reaction_mat(self):
        if not (self.sx_mat.shape == (self.num_joints, 3) and self.sy_mat.shape == (self.num_joints,3)):
            print("Invalid reaction matrices!!")
            exit(127)
    
    def check_valid_l_mat(self):
        if not (self.l_vect.shape == (self.num_joints*2,1)):
            print("Invalid load vector!!")
            exit(127)
    
    def check_valid_pos_vect(self):
        if not (self.x_vect.shape == (1,self.num_joints) and self.y_vect.shape == (1,self.num_joints)):
            print("Invalid joint position vectors!!")
            exit(127)
    
                
##### Mutators #####
    def set_truss_design(self, c, sx, sy, x, y, l):
        self.c_mat = c
        self.sx_mat = sx
        self.sy_mat = sy
        self.x_vect = x
        self.y_vect = y
        self.l_vect = l
        self.load = abs(np.sum(self.l_vect))
        self.num_joints = x.size
        self.num_mems = self.num_joints * 2 - 3
        self.check_valid_truss()
        self.check_valid_reaction_mat()
        self.check_valid_pos_vect()
        self.check_valid_l_mat()
        self.check_valid_c_mat()

        self.len_vect = []
        self.p_crit_vect = []
        for i in range(self.num_mems):
            j1_pos, j2_pos = self.find_joint_pos_for_member(i)
            self.len_vect.append(self.calc_joint_dist(j1_pos, j2_pos))
        self.p_crit_vect = [self.calc_buckling_f(l) for l in self.len_vect]
        self.total_mem_len = sum(self.len_vect)
        self.weight_vect = [self.weight_per_oz * l for l in self.len_vect]

    def set_c_mat(self, c):
        self.c_mat = c

    def set_sx_mat(self, sx):
        self.sx_mat = sx
    
    def set_sy_mat(self, sy):
        self.sy_mat = sy
    
    def set_x_vect(self, x):
        self.x_vect = x
    
    def set_y_vect(self, y):
        self.y_vect = y
    
    def set_l_vect(self, l):
        self.l_vect = l
    
    def set_num_joints(self, joints):
        self.num_joints = joints

    def set_num_mems(self, members):
        self.num_mems = members
    
    def set_total_mem_len(self, mem_len):
        self.total_mem_len = mem_len
    
##### Accessors #####

    def get_c_mat(self):
        return self.c_mat

    def get_sx_mat(self):
        return self.sx_mat
    
    def get_sy_mat(self):
        return self.sy_mat
    
    def get_x_vect(self):
        return self.x_vect
    
    def get_y_vect(self):
        return self.y_vect
    
    def get_l_vect(self):
        return self.l_vect
    
    def get_num_joints(self):
        return self.num_joints 

    def get_num_mems(self):
        return self.num_mems
    
    def get_total_mem_len(self):
        return self.total_mem_len
    
    def get_a_mat(self):
        return self.a_mat
    

