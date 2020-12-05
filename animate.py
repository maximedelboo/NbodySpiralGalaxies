import matplotlib; matplotlib.use("TkAgg")
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from matplotlib import style
import numpy as np
from multiprocessing import Process, Pipe, Queue
import model_test


def run_animation(queue):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    ani = animation.FuncAnimation(fig, animate, interval=100, fargs=(ax, queue,))
    plt.show()


def animate(i, ax, queue):
    # r_x = np.array(open('model_x.txt', 'r').read().split('\n'))[:-1].astype(np.float)
    # r_y = np.array(open('model_y.txt', 'r').read().split('\n'))[:-1].astype(np.float)

    if not queue.empty():
        r_x, r_y = queue.get()
        ax.clear()
        ax.plot(r_x, r_y, ls='', marker='.')
        ax.set_aspect('equal')

if __name__ == '__main__':
    queue = Queue()

    model_process = Process(target=model_test.model_loop, args=(queue,))
    animate_process = Process(target=run_animation, args=(queue,))

    model_process.start()
    animate_process.start()