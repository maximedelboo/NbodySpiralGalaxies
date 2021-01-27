from multiprocessing.sharedctypes import RawArray

import pandas as pd
import numpy as np
import multiprocessing as mp
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time

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
        self.data = None
        self.clearer = None

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


def read_data(queue, l_queue):
    filename = 'C:/Users/tip/Files/Git projects/Planeten-enal/huts d9.csv'
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


def wait_for_load(l_queue, f_queue):
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
    f_queue.put([data])
    rd2.data = data
    print("Finished loading!")


def clear_queue(queue):
    while not queue.empty():
        try:
            queue.get(False)
        except queue.Empty:
            continue
        queue.task_done()
    queue.join()


def run_animation(queue, f_queue):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set(xlim=(-1e18, 1e18), ylim=(-1e18, 1e18))

    while queue.qsize() < 3:
        pass
    reading_chunk = queue.get(block=True)
    row1 = reading_chunk[0]
    row1 = row1.reshape(int(len(row1) / 2), 2).transpose()

    line = ax.plot(row1[0], row1[1], '.', ms=1)[0]
    rd = ReadData(reading_chunk)
    rd2.set_reading_chunk(reading_chunk)

    queue.task_done()

    ani = animation.FuncAnimation(fig, animate, interval=30, frames=iterations, fargs=(ax, line, rd, queue, f_queue))
    plt.draw()
    plt.show()


def animate(i, ax, line, rd, queue, f_queue):
    if i < iterations - 1:
        frame_i = i + 1
    else:
        frame_i = i % iterations

    ax.set_title(str(frame_i))
    start_time = time.perf_counter()
    if not rd.finished_loading and not f_queue.empty():
        rd.data = f_queue.get(block=True)[0]
        time0 = time.perf_counter() - start_time
        start_time = time.perf_counter()
        time1 = time.perf_counter() - start_time
        start_time = time.perf_counter()
        rd.finished_loading = True
        f_queue.task_done()
        f_queue.join()
        time2 = time.perf_counter() - start_time
        print("Put data!")
        clear_process = mp.Process(target=clear_queue, args=(queue,))
        clear_process.start()
        rd.clearer = clear_process
        print(f"{time0:0.8f}" + " " + f"{time1:0.8f}" + " " + f"{time2:0.8f}")

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

    line.set_xdata(row[0])
    line.set_ydata(row[1])


def init_animation(row1):
    poii_num = 4
    for i in range(poii_num):
        poii = np.random.randint(int(len(row1) / 2))
    poii = np.random.randint(int(len(row1) / 2))
    poii2 = np.random.randint(int(len(row1) / 2))


if __name__ == '__main__':
    main_queue = mp.Manager().Queue()
    loaded_queue = mp.Manager().Queue()
    finished_queue = mp.Manager().Queue()

    processes = []

    rd2 = ReadData()

    read_process = mp.Process(target=read_data, args=(main_queue, loaded_queue,))
    processes.append(read_process)
    wait_load_process = mp.Process(target=wait_for_load, args=(loaded_queue, finished_queue,))
    processes.append(wait_load_process)
    animate_process = mp.Process(target=run_animation, args=(main_queue, finished_queue,))
    processes.append(animate_process)

    read_process.start()
    animate_process.start()
    wait_load_process.start()

    for i, process in enumerate(processes):
        process.join()
