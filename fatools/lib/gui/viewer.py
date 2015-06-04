
import tkinter

from fatools.lib.fautil.wavelen2rgb import wavelen2rgb

class SimplePlot(tkinter.Canvas):

    def plot_channel(self, channel):
        
        gram = channel.raw
        color = '#%02x%02x%02x' % tuple(wavelen2rgb(channel.wavelength))

        y0 = gram[0]

        for i in range(1, len(gram)):
            y1 = gram[i]
            self.create_line( i-1, y0, i, y1, fill=color )
            y0 = y1
        



def viewer( trace ):

    root = tkinter.Tk()

    widget = SimplePlot(root)
    widget.pack(fill='both', expand=1)
    widget.config(background='#000000')

    channels = trace.get_channels()

    for dye in channels:
        widget.plot_channel( channels[dye] )

    widget.update()

    root.mainloop()


def viewer( trace ):

    from matplotlib import pyplot as plt

    channels = trace.get_channels()

    for dye in channels:
        c = channels[dye]
        plt.plot( c.smooth(), color= tuple( x/255 for x in wavelen2rgb(c.wavelength)) )

    plt.show()
