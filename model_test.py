import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.animation as animation
import scipy.stats as sps

r_x = list()
r_y = list()
n_radius = 50
rough_n = 50000

max_radius = 25e19

radii = np.arange(0, 50, 1)
radii = max_radius / n_radius * (radii + 1)
circumferences = 2 * np.pi * radii
total_circumference = np.sum(circumferences)

for i in range(n_radius):
    circumference = circumferences[i]
    radius = radii[i]
    # j_tot = int(np.floor(total_circumference / rough_n * circumference))
    j_tot = int(np.floor(rough_n / (n_radius**2) * (np.random.random() + 0.5)))
    for j in range(j_tot):
        phi = np.random.random() * 2 * np.pi
        r_x.append(np.cos(phi) * radius)
        r_y.append(np.sin(phi) * radius)

n = len(r_x)
r_x = np.array(r_x)
r_y = np.array(r_y)
rd_x = np.zeros(n)
rd_y = np.zeros(n)
rdd_x = np.zeros(n)
rdd_y = np.zeros(n)
m = np.random.gamma(shape=0.7, scale=0.47, size=n) + 0.05

cells = list()
theta = 0.98



def M(r, alpha, r_0, rho_0):

    return 2*np.pi*rho_0*(r_0**alpha)*(r**(2-alpha))/(2-alpha)

#hoeksnelheid door darkmatter halo
def dphi(r, alpha, r_0, omega_0):
    
    return omega_0*((r/r_0)**(-alpha/2))

class Cell:
    def __init__(self, top_left_x, top_left_y, size, parent=None):
        self.x = top_left_x
        self.y = top_left_y
        self.size = size
        self.parent = parent
        self.children = None
        self.point_x = None
        self.point_y = None
        self.point_m = None
        self.cm_x = None
        self.cm_y = None
        self.total_m = None

    def quadruple(self):
        global cells
        new_size = self.size / 2
        self.children = list()
        x_new_tl = self.x
        y_new_tl = self.y
        x_new_tr = x_new_tl + new_size
        y_new_tr = y_new_tl
        x_new_bl = x_new_tl
        y_new_bl = y_new_tl - new_size
        x_new_br = x_new_tl + new_size
        y_new_br = y_new_tl - new_size
        cell_tl = Cell(x_new_tl, y_new_tl, new_size, self)
        cell_tr = Cell(x_new_tr, y_new_tr, new_size, self)
        cell_bl = Cell(x_new_bl, y_new_bl, new_size, self)
        cell_br = Cell(x_new_br, y_new_br, new_size, self)
        self.children = [cell_tl, cell_tr, cell_bl, cell_br]
        cells += [cell_tl, cell_tr, cell_bl, cell_br]

    def find_child_cells(self, new_point_x, new_point_y, new_point_m):
        for child in self.children:
            if self.point_x and child.x <= self.point_x <= child.x + child.size and \
                    child.y - child.size <= self.point_y <= child.y:
                child.insert(self.point_x, self.point_y, self.point_m)
                self.point_x = None
                self.point_y = None
                self.point_m = None
            if new_point_x and child.x <= new_point_x <= child.x + child.size and \
                    child.y - child.size <= new_point_y <= child.y:
                child.insert(new_point_x, new_point_y, new_point_m)
                new_point_x = None

    def insert(self, new_point_x, new_point_y, new_point_m):
        if self.point_x:
            self.quadruple()
            self.find_child_cells(new_point_x, new_point_y, new_point_m)
        elif self.children:
            self.find_child_cells(new_point_x, new_point_y, new_point_m)
        else:
            self.point_x = new_point_x
            self.point_y = new_point_y
            self.point_m = new_point_m

    def compute_mass(self):
        if self.point_x:
            self.cm_x = self.point_x
            self.cm_y = self.point_y
            self.total_m = self.point_m
        elif self.children:
            children = self.children
            c1_cm_x, c1_cm_y, c1_m = children[0].compute_mass()
            c2_cm_x, c2_cm_y, c2_m = children[1].compute_mass()
            c3_cm_x, c3_cm_y, c3_m = children[2].compute_mass()
            c4_cm_x, c4_cm_y, c4_m = children[3].compute_mass()
            total_m = (c1_m + c2_m + c3_m + c4_m)
            self.cm_x = (c1_cm_x * c1_m + c2_cm_x * c2_m + c3_cm_x * c3_m + c4_cm_x * c4_m) / total_m
            self.cm_y = (c1_cm_y * c1_m + c2_cm_y * c2_m + c3_cm_y * c3_m + c4_cm_y * c4_m) / total_m
            self.total_m = total_m
        else:
            return 0, 0, 0
        return self.cm_x, self.cm_y, self.total_m

    def get_mass_distance_term(self, probe_x, probe_y):
        if self.point_x:
            x = self.point_x
            y = self.point_y
            d = np.sqrt((x - probe_x) ** 2 + (y - probe_y) ** 2)
            if d == 0:
                return 0, 0
            md_x = self.point_m * (x - probe_x) / np.power(d, 3)
            md_y = self.point_m * (y - probe_y) / np.power(d, 3)

            return md_x, md_y
        elif self.children:
            global theta
            x = self.cm_x
            y = self.cm_y
            d = np.sqrt((x - probe_x) ** 2 + (y - probe_y) ** 2)
            if d == 0:
                return 0, 0
            if self.size / d < theta:
                md_x = self.total_m * (x - probe_x) / np.power(d, 3)
                md_y = self.total_m * (y - probe_y) / np.power(d, 3)
                return md_x, md_y
            else:
                md_x = 0
                md_y = 0
                for child in self.children:
                    p_md_x, p_md_y = child.get_mass_distance_term(probe_x, probe_y)
                    md_x += p_md_x
                    md_y += p_md_y
                return md_x, md_y
        else:
            return 0, 0

    def draw(self):
        # ax.add_patch(Rectangle((self.x, self.y - self.size), self.size, self.size, fill=None))
        pass


G = 6.67408e-11
m_sol = 1.98847e30
dt = 31e12

# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1)

quad_tree = Cell(-max_radius * 5, max_radius * 5, 10 * max_radius)

for i in range(n):
    quad_tree.insert(r_x[i], r_y[i], m[i])

quad_tree.compute_mass()

for i in range(n):
    r_x_i_0 = r_x[i]
    r_y_i_0 = r_y[i]
    l_i_factor = ((2 * np.random.random() + 1) * 1e5) / np.sqrt(r_x_i_0 ** 2 + r_y_i_0 ** 2) * 0.1
    rd_x_i_0 = r_y_i_0 * l_i_factor
    rd_y_i_0 = -r_x_i_0 * l_i_factor
    mdt_x_i_0, mdt_y_i_0 = quad_tree.get_mass_distance_term(r_x_i_0, r_y_i_0)
    rdd_x_i_0 = -G * m_sol * mdt_x_i_0
    rdd_y_i_0 = -G * m_sol * mdt_y_i_0
    r_x_i_1 = r_x_i_0 + rd_x_i_0 * dt + 0.5 * rdd_x_i_0 * dt ** 2
    r_y_i_1 = r_y_i_0 + rd_x_i_0 * dt + 0.5 * rdd_y_i_0 * dt ** 2
    r_x[i] = r_x_i_1
    r_y[i] = r_y_i_1
    rd_x[i] = rd_x_i_0
    rd_y[i] = rd_y_i_0
    rdd_x[i] = rdd_x_i_0
    rdd_y[i] = rdd_y_i_0


# Notice that r is 1 ahead of rd and rdd


def model_loop(queue):
    #dark matter halo
    alpha = 1.5
    m_sol = 1.98847e30
    f= 0.2
    M_d = m_sol*sum(m)
    M_h = M_d/f - M_d
    R = max_radius
    r_0 = max_radius
    rho_0 = (2-alpha)*M_h/(2*np.pi*(r_0**alpha)*(R**(2-alpha)))
    
    G = 6.67408e-11
    omega_0 = np.sqrt(2*np.pi*G*rho_0/(2-alpha))
    
    for j, t in enumerate(np.arange(0, dt * 10000, dt), 1):
        print(t)

        cells = list()
        # plt.plot(r_x, r_y, ls='', marker='.')
        quad_tree = Cell(-max_radius * 10, max_radius * 10, 20 * max_radius)

        for i in range(n):
            quad_tree.insert(r_x[i], r_y[i], m[i])

        # for cell in cells:
        #    cell.draw()

        quad_tree.compute_mass()

        for i in range(n):
            r_x_i_j = r_x[i]
            r_y_i_j = r_y[i]
            rd_x_i_jm1 = rd_x[i]
            rd_y_i_jm1 = rd_y[i]
            rdd_x_i_jm1 = rdd_x[i]
            rdd_y_i_jm1 = rdd_y[i]
            
            r = np.sqrt(r_x_i_j**2 + r_y_i_j**2)
            

            
            a_x = (r_x_i_j/r)*G*M(r, alpha, r_0, rho_0)/r**2
            a_y = (r_y_i_j/r)*G*M(r, alpha, r_0, rho_0)/r**2

            # calculate rdd_i_j
            md_x_i_j, md_y_i_j = quad_tree.get_mass_distance_term(r_x_i_j, r_y_i_j)
            rdd_x_i_j = -G * m_sol * md_x_i_j + a_x
            rdd_y_i_j = -G * m_sol * md_y_i_j + a_y

            # calculate rd_i_j
            rd_x_i_j = rd_x_i_jm1 + 0.5 * (rdd_x_i_jm1 + rdd_x_i_j) * dt
            rd_y_i_j = rd_y_i_jm1 + 0.5 * (rdd_y_i_jm1 + rdd_y_i_j) * dt

            # calculate r_i_(j+1)
            r_x_i_jp1 = r_x_i_j + rd_x_i_j * dt + 0.5 * rdd_x_i_j * dt ** 2
            r_y_i_jp1 = r_y_i_j + rd_y_i_j * dt + 0.5 * rdd_y_i_j * dt ** 2

            r_x[i] = r_x_i_jp1
            r_y[i] = r_y_i_jp1
            rd_x[i] = rd_x_i_j
            rd_y[i] = rd_y_i_j
            rdd_x[i] = rdd_x_i_j
            rdd_y[i] = rdd_y_i_j

        # np.savetxt('model_x.txt', r_x)
        # np.savetxt('model_y.txt', r_y)

        queue.put((r_x, r_y))

        # ax.set_aspect('equal')
        # plt.show()
