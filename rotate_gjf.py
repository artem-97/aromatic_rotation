import numpy as np
from scipy.spatial.transform import Rotation as R

# script for Gaussian09 input files generation

table = {
    1: 'H',
    6: 'C',
    7: 'N',
    8: 'O',
    16: 'S',
}

filename = 'furan.xyz'
n_template = 46  # atoms in template
n_points = 12  # number of rotattion points
header = '''%mem=2GB
#BP86/Def2SVP Opt EmpiricalDispersion(GD3BJ) scrf=(cpcm,solvent=acetonitrile) Freq

furan

0 1'''

with open(filename) as file:
    template = []
    substrate = []
    cnt = n_template
    for line in file:
        if cnt > 0:
            cnt -= 1
            splitted = line.split()
            charge = int(splitted[0])
            x = float(splitted[1])
            y = float(splitted[2])
            z = float(splitted[3])
            template.append([charge, np.array([x, y, z])])
        else:
            splitted = line.split()
            charge = int(splitted[0])
            x = float(splitted[1])
            y = float(splitted[2])
            z = float(splitted[3])
            substrate.append([charge, np.array([x, y, z])])
    Rc = np.array([0., 0., 0.])
    M = 0
    for a in substrate:
        M += 1
        Rc += a[1]
    Rc /= M
    for i in range(len(substrate)):
        substrate[i][1] -= Rc
    for i in range(len(template)):
        template[i][1] -= Rc
    A = substrate[0][1]
    B = substrate[1][1]
    C = substrate[2][1]
    u = A - C
    v = B - C
    n = np.cross(u, v)
    n /= np.sqrt(np.sum(n**2))
    for k in range(n_points):
        phi = k * 2 * np.pi / n_points
        r = R.from_rotvec(phi * n).as_matrix()
        for i in range(len(substrate)):
            substrate[i][1] = r @ substrate[i][1]
        file_k = str(k) + '.gjf'
        with open(file_k, 'w') as res:
            print(header, file=res)
            for e in template:
                print(table[e[0]], e[1][0], e[1][1], e[1][2], file=res)
            for e in substrate:
                print(table[e[0]], e[1][0], e[1][1], e[1][2], file=res)
            print("", file=res)
