from Analysis import Analysis
import numpy as np

def main():
    A = Analysis(8, 13, 20)
    c = np.array([[ 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [  1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                 [  0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0],
                 [  0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
                 [  0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0],
                 [  0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0],
                 [  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1],
                 [  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]])
    sx = np.array([ [1, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0]])
    sy = np.array([ [0, 1, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 1]])
    x = np.array([0, 0, 4, 4, 8, 8, 12, 12])
    y = np.array([0, 4, 4, 8, 8, 4, 4, 0])
    x = x.reshape(1,8)
    y = y.reshape(1,8)

    l = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -25, 0, 0, 0, 0, 0])
    l = l.reshape(16,1)
    A.set_truss_design(c, sx, sy, x, y, l, 8, 13, 20)
    A.mat_save()
    A.mat_load()
    A.construct_A_mat()
    # A.solve()
    A.print_output()
    # A.print_truss_details()

    # print(A.calc_joint_dist([3,0], [0,0]))
    
    
if __name__ == "__main__":
    main()