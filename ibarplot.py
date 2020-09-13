class InteractiveBarPlot:
    
    def __init__(self, title='Init title', colormap='Spectral', colorslice=[0.09, 0.89], yrange=[0, 10], xrange=[0, 10]):
        
        fig, ax = plt.subplots()
        self.figure = fig
        self.ax = ax
        
        self.xrange = xrange
        self.yrange = yrange
        
        self.ax.set_xlim(xrange)
        self.ax.set_ylim(yrange)
        self.ax.set_title(title)
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        
        self.__FACTOR = 0.95
        
        self.colormap = cm.get_cmap(name=colormap)
        self.mincolor = colorslice[0]
        self.maxcolor = colorslice[1]
        self.coloredges = (self.colormap(self.mincolor), self.colormap(self.mincolor), self.colormap(self.maxcolor))
        
        self.values = []
        self.colors = []
        self.err = []
        self.edgecolors = []
        self.widths = []
        self.widths_scale = []
        self.positions = []
        self.indexes = []
        
        self.__delta = 0
        self.mean = 0
        self.top_lim = (yrange[1] - yrange[0])*0.9
        self.bot_lim = (yrange[1] - yrange[0])*0.1
        self.common_width = 0
        self.bar_graph = None
        self.bar_graph_bars = None
        self.line_top = None
        self.line_bot = None
        self.fill_between = None
        self.scat_top = None
        self.scat_bot = None
        self.graph_labels = None
        self.labels_data = None
        self.top_clicked = False
        self.bot_clicked = False
        self.cent_clicked = False
        
        
        self.new_width = lambda l: self.common_width*l/np.sum(self.widths_scale)
        self.new_pos = lambda i: np.sum(self.widths[:i]) + np.mean([self.widths[i], self.widths[0]]) + (i + 1)*0.2
        
        self.vect_update_colors = np.vectorize(self.__update_color, otypes=[list])
        self.vect_update_width = np.vectorize(self.new_width, otypes=[list])
        self.vect_update_pos = np.vectorize(self.new_pos, otypes=[list])
        self.vect_add_labels = np.vectorize(self.__add_bar_label, otypes=[list])
    
    def set_lims(self, top_lim=-1, bot_lim=1):
        
        self.top_lim = top_lim
        self.bot_lim = bot_lim
        self.mean = np.mean([top_lim, bot_lim])
        
    def __change_color(self, color, amount=0.5):
    
        try:
            c = mc.cnames[color]
        except:
            c = color
            
        c = rgb_to_hls(*mc.to_rgb(c))
        
        return hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])
        
    def add_data(self, values=[], err=[], autoscale=True, labels=None):
        
        self.values = np.array(values)
        self.err = np.array(err)
        self.__delta = np.max(self.values) * 0.02
        size = self.values.size
        
        self.colors = [0 for i in range(size)]
        self.edgecolors = [0 for i in range(size)]
        self.widths = [0 for i in range(size)]
        self.widths_scale = [0 for i in range(size)]
        self.indexes = np.array([i for i in range(size)])
        self.common_width = size
        
        if autoscale:
            
            self.xrange = [0, round(size*1.5, 0)]
            self.ax.set_xlim(self.xrange)
            
        if not labels is None and np.array(labels).size == size:
            
            self.graph_labels = labels
            self.labels_data = [0 for i in range(size)]
        
    def __update_color(self, index):
        
        scale = 1
        top = max(self.bot_lim, self.top_lim)
        bot = min(self.bot_lim, self.top_lim)
        self.mean = np.mean([top, bot])
        mean = self.mean
        
        value = self.values[index]
        
        if bot < value < top:
            
            dist = (value - bot) / (top - bot)
            scale += dist if dist <= 0.5 else 1 - dist
            color = self.colormap(self.mincolor + dist*(self.maxcolor - self.mincolor))

        else:
            
            sign = np.sign(bot - value).astype('int64')
            color = self.coloredges[sign]
            
        self.colors[index] = color
        self.edgecolors[index] = self.__change_color(color, 0.7)
        self.widths_scale[index] = scale
    
    def __update_colors(self):
        
        self.vect_update_colors(self.indexes)
        
        
        self.widths = self.vect_update_width(self.widths_scale)
        self.positions = self.vect_update_pos(self.indexes)
    
    def __add_bar_label(self, index):
        
        bar = self.bar_graph_bars[index]
        height = bar.get_height()
        text = self.graph_labels[index]
        ratio = 1 if len(text) == 1 else 1/np.log(len(text))
        fontsize = self.widths[index]*15*ratio + 0.4
        
        self.labels_data[index] = self.ax.text(bar.get_x() + bar.get_width()/2.0, 0.5*height,
                                              text, ha='center', va='bottom',
                                              color=self.__change_color(self.colors[index], 0.6),
                                              fontsize=fontsize)
    
    def __add_bar_labels(self):
        
        if not self.graph_labels is None:
            
            self.vect_add_labels(self.indexes)
            
            
        
    def __plot(self, draw_tline=True, draw_bline=True):
        
        self.__update_colors()
        
        self.bar_graph = self.ax.bar(self.positions, self.values, yerr=self.err,
                                    capsize=16, color=self.colors, linewidth=1,
                                    edgecolor=self.edgecolors, width=self.widths)
        
        indexed_bars = [bar for i, bar in enumerate(self.bar_graph)]
        self.bar_graph_bars = indexed_bars
        
        self.__add_bar_labels()
        
        self.fill_between = self.ax.fill_between([0, self.xrange[1]*0.9], [self.top_lim, self.top_lim], 
                                                 [self.bot_lim, self.bot_lim], color='skyblue', alpha=0.15)
        
        if draw_tline:
            
            self.line_top, = self.ax.plot([0, self.xrange[1]*0.9], [self.top_lim, self.top_lim], 
                                          color='skyblue', linewidth=1.5, alpha=0.4)
            
            self.scat_top = self.ax.scatter([self.xrange[1]*0.9], [self.top_lim], marker='o', 
                                            color='skyblue', s=80, zorder=10)
        
        if draw_bline:
            
            self.line_bot, = self.ax.plot([0, self.xrange[1]*0.9], [self.bot_lim, self.bot_lim], 
                                          color='skyblue', linewidth=1.5, alpha=0.4)
            
            self.scat_bot = self.ax.scatter([self.xrange[1]*0.9], [self.bot_lim], marker='o', 
                                            color='skyblue', s=80, zorder=10)
    
    def __connect_events(self):
        
        self.figure.canvas.mpl_connect('button_press_event', self.__on_drag_start)
        self.figure.canvas.mpl_connect('motion_notify_event', self.__drag_event)
        self.figure.canvas.mpl_connect('button_release_event', self.__on_drag_end)
        self.figure.canvas.mpl_connect('axes_leave_event', self.__on_drag_end)
    
    def __on_drag_start(self, event):
        
        shift = self.__FACTOR * (self.mean - self.top_lim)
        
        if self.mean - shift < event.ydata < self.mean + shift:
            self.cent_clicked = True
        
        elif self.top_lim - self.__delta < event.ydata < self.top_lim + self.__delta:
            self.top_clicked = True
        
        elif self.bot_lim - self.__delta < event.ydata < self.bot_lim + self.__delta:
            self.bot_clicked = True
            
    
    def __drag_event(self, event):
        
        if self.top_clicked:
            
            self.line_top.remove()
            self.fill_between.remove()
            self.scat_top.remove()
            self.bar_graph.errorbar.remove()
            _ = [p.remove() for p in self.bar_graph.patches]
            
            if self.labels_data: _ = [l.remove() for l in self.labels_data]
            
            self.top_lim = event.ydata
            
            self.__plot(draw_bline=False)
        
        if self.bot_clicked:
            
            self.line_bot.remove()
            self.fill_between.remove()
            self.scat_bot.remove()
            self.bar_graph.errorbar.remove()
            _ = [p.remove() for p in self.bar_graph.patches]
            if self.labels_data: _ = [l.remove() for l in self.labels_data]
                
            self.bot_lim = event.ydata
            
            self.__plot(draw_tline=False)
        
#         if self.cent_clicked:
            
#             self.line_bot.remove()
#             self.line_top.remove()
#             self.fill_between.remove()
#             self.scat_bot.remove()
#             self.scat_top.remove()
#             self.bar_graph.errorbar.remove()
#             _ = [p.remove() for p in self.bar_graph.patches]
#             if self.labels_data: _ = [l.remove() for l in self.labels_data]

            
            
    def __on_drag_end(self, event):
        
        self.cent_clicked = False
        self.top_clicked = False
        self.bot_clicked = False
        
    
    def launch(self):
        
        self.__connect_events()
        self.__plot()