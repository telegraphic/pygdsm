def show_plt(): # pragma: no cover
    try:
        plt.show() # noqa: F821
    except:
        import pylab as plt
        plt.show()