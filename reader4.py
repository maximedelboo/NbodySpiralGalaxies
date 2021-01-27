import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from multiprocessing import Process, Manager
from multiprocessing.shared_memory import SharedMemory

matplotlib.use("TkAgg")

chunksize = 100
iterations = 3000
fps = 60


class ReadData:
    def __init__(self, reading_chunk=None):
        self.reading_chunk = reading_chunk
        self.is_loading = False
        self.finished_loading = False
        self.data_array = None
        self.data_shared = None
        self.data = None
        self.clearer = None
        self.poiis = None
        self.poiis_x = None
        self.poiis_y = None

    def set_reading_chunk(self, reading_chunk):
        self.reading_chunk = reading_chunk

    def get_row(self, i):
        if not self.finished_loading:
            return self.reading_chunk[i % chunksize]
        else:
            if self.clearer and not self.clearer.is_alive():
                self.clearer.join()
                self.clearer = None
                print("Cleared!")
            return self.data[i]

    def update_poii_data(self, transposed_row):
        for i, x_list in enumerate(self.poiis_x):
            # print("appendi " + str(transposed_row[0][self.poiis[i]]))
            x_list.append(transposed_row[0][self.poiis[i]])
            self.poiis_y[i].append(transposed_row[1][self.poiis[i]])



def read_data(queue, l_queue):
    filename = 'C:/Users/tip/Files/Git projects/Planeten-enal/huts d6p5 no pp 2p9r.csv'
    # path van de csv, best bij de python file
    reader = pd.read_csv(filename, sep=";", chunksize=chunksize, header=None)
    i = 0
    for chunk in reader:
        print(i * chunksize)
        data_chunk = chunk.to_numpy()[:, :-1]
        queue.put(data_chunk)
        l_queue.put(data_chunk)
        i += 1
    num_chunks = i
    l_queue.put(num_chunks)


def loaded_data(l_queue, f_queue):
    rows = []

    finished_loading = False
    while not finished_loading:
        next_chunk = l_queue.get(block=True)
        if isinstance(next_chunk, int):
            finished_loading = True
        else:
            for row in next_chunk:
                rows.append(row)
        l_queue.task_done()

    data = np.array(rows)
    l_queue.join()
    f_queue.put(data)
    print("Finished loading!")


def clear_queue(queue):
    while not queue.empty():
        try:
            queue.get(False)
        except queue.Empty:
            continue
        queue.task_done()
    queue.join()


def run_animation(queue, a_queue):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set(xlim=(-1e18, 1e18), ylim=(-1e18, 1e18))

    while queue.qsize() < 2:
        pass
    reading_chunk = queue.get(block=True)
    row1 = reading_chunk[0]
    row1 = row1.reshape(int(len(row1) / 2), 2).transpose()

    line = ax.plot(row1[0], row1[1], '.', ms=1)[0]
    rd = ReadData(reading_chunk)

    queue.task_done()

    Gvar = 6.67408e-11
    m_sol = 1.98847e30
    m_center_mass = 700
    radius_gal = 1e18
    omega = np.sqrt(Gvar * m_center_mass * m_sol / (radius_gal / 2) ** 3)
    dt = 35e12
    p_c = np.linspace(0, 1, 100) * 2 * np.pi
    r_in = (Gvar * m_center_mass * m_sol / (4 * omega ** 2)) ** (1 / 3)
    r_out = (9 * Gvar * m_center_mass * m_sol / (4 * omega ** 2)) ** (1 / 3)
    x_c_i = np.cos(p_c) * r_in
    y_c_i = np.sin(p_c) * r_in
    x_c_o = np.cos(p_c) * r_out
    y_c_o = np.sin(p_c) * r_out

    constants = {'Gvar': Gvar, 'm_sol': m_sol, 'm_center_mass': m_center_mass, 'radius_gal': radius_gal, 'omega': omega,
                 'dt': dt}
    lines_poiis = get_lines_poiis(row1, ax, radius_gal)
    rd.poiis = lines_poiis['poiis']
    rd.poiis_x = [[] for i in range(len(lines_poiis['poiis']))]
    rd.poiis_y = [[] for i in range(len(lines_poiis['poiis']))]

    res_line_i = ax.plot(x_c_i, y_c_i, ms=5, c='purple')[0]
    res_line_o = ax.plot(x_c_o, y_c_o, ms=5, c='purple')[0]

    ani = animation.FuncAnimation(fig, animate, interval=30, frames=iterations, fargs=(ax, line, rd, queue, a_queue,
                                                                                       lines_poiis, constants))
    plt.draw()
    plt.show()


def animate(i, ax, line, rd, queue, a_queue, lines_poiis, constants):
    if i < iterations - 1:
        frame_i = i + 1
    else:
        frame_i = i % iterations

    ax.set_title(str(frame_i))
    start_time = time.perf_counter()

    if not rd.finished_loading and not a_queue.empty():
        data_shared_info = a_queue.get(block=True)
        rd.data_shared = SharedMemory(data_shared_info['name'])
        rd.data = np.ndarray(shape=data_shared_info['shape'], dtype=data_shared_info['dtype'], buffer=rd.data_shared.buf)
        time0 = time.perf_counter() - start_time
        start_time = time.perf_counter()
        time1 = time.perf_counter() - start_time
        start_time = time.perf_counter()
        rd.finished_loading = True
        a_queue.task_done()
        a_queue.join()
        time2 = time.perf_counter() - start_time
        print("Put data!")
        # clear_process = Process(target=clear_queue, args=(queue,))
        # clear_process.start()
        # rd.clearer = clear_process
        # print(f"{time0:0.8f}" + " " + f"{time1:0.8f}" + " " + f"{time2:0.8f}")


    start_time = time.perf_counter()
    if not rd.finished_loading and frame_i % chunksize == 0:
        rd.set_reading_chunk(queue.get(block=True))
        queue.task_done()
    time4 = time.perf_counter() - start_time
    start_time = time.perf_counter()
    row = rd.get_row(frame_i)
    row = row.reshape(int(len(row) / 2), 2).transpose()
    time5 = time.perf_counter() - start_time
    # print(f"{time4:0.8f}" + " " + f"{time5:0.8f}")

    t = frame_i*constants['dt']
    T = constants['omega']*t
    o_x = constants['radius_gal'] * np.cos(T)
    o_y = constants['radius_gal'] * np.sin(T)
    lines_poiis['o_line'].set_xdata([0, o_x])
    lines_poiis['o_line'].set_ydata([0, o_y])

    rd.update_poii_data(row)
    x_poiis = rd.poiis_x
    y_poiis = rd.poiis_y
    s_i = min(len(rd.poiis_x[0]), 300)
    for i, poii_line in enumerate(lines_poiis['lines']):
        poii_line.set_xdata(x_poiis[i][-s_i:])
        poii_line.set_ydata(y_poiis[i][-s_i:])

    line.set_xdata(row[0])
    line.set_ydata(row[1])


def get_lines_poiis(row1, ax, radius_gal):
    poii_num = 0
    poiis = []
    for i in range(poii_num):
        poii = np.random.randint(int(len(row1[0]) / 2))
        poiis.append(poii)
    extra_poiis = [4510]
    all_poiis = poiis + extra_poiis
    lines = []
    for poii in all_poiis:
        line = ax.plot(row1[0][poii], row1[1][poii], ms=10)[0]
        lines.append(line)

    omega_line = ax.plot([0, radius_gal], [0, 0], ms=5, c='black')[0]

    return {'lines': lines, 'poiis': all_poiis, 'o_line': omega_line}

if __name__ == '__main__':
    main_queue = Manager().Queue()
    loaded_queue = Manager().Queue()
    finished_queue = Manager().Queue()
    address_queue = Manager().Queue()

    processes = []

    read_process = Process(target=read_data, args=(main_queue, loaded_queue,))
    processes.append(read_process)
    animate_process = Process(target=run_animation, args=(main_queue, address_queue,))
    processes.append(animate_process)
    loaded_data_process = Process(target=loaded_data, args=(loaded_queue, finished_queue,))
    processes.append(loaded_data_process)

    read_process.start()
    animate_process.start()
    loaded_data_process.start()

    finished_data = finished_queue.get(block=True)
    shp_data = finished_data.shape
    dtype_data = finished_data.dtype
    data_shm = SharedMemory(create=True, size=finished_data.nbytes)
    finished_data_shared = np.ndarray(shape=shp_data, dtype=dtype_data, buffer=data_shm.buf)
    finished_data_shared[:] = finished_data[:]
    finished_queue.task_done()
    finished_queue.join()

    address_queue.put({'name': data_shm.name, 'shape': shp_data, 'dtype': dtype_data})

    for p_i, process in enumerate(processes):
        process.join()
