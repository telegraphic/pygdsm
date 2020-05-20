def show_plt():
    try:
        plt.show()
    except:
        import pylab as plt
        plt.show()