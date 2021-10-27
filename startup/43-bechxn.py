from bluesky.callbacks.best_effort import BestEffortCallback
from bluesky.callbacks.best_effort import LivePlotPlusPeaks
import hxntools.scans

from bluesky.callbacks.core import CallbackBase, LiveTable
from bluesky.callbacks.mpl_plotting import LivePlot, LiveGrid, LiveScatter

import matplotlib.pyplot as plt

from itertools import count

def hinted_fields(descriptor):
    # Figure out which columns to put in the table.
    obj_names = list(descriptor['object_keys'])
    # We will see if these objects hint at whether
    # a subset of their data keys ('fields') are interesting. If they
    # did, we'll use those. If these didn't, we know that the RunEngine
    # *always* records their complete list of fields, so we can use
    # them all unselectively.
    columns = []
    for obj_name in obj_names:
        try:
            fields = descriptor.get('hints', {}).get(obj_name, {})['fields']
        except KeyError:
            fields = descriptor['object_keys'][obj_name]
        columns.extend(fields)
    return columns

class LiveScatterHxn(LiveScatter):

    def start(self, doc):
        self._xdata.clear()
        self._ydata.clear()
        self._Idata.clear()
        sc = self.ax.scatter(self._xdata, self._ydata, 
                             norm=self._norm, cmap=self.cmap, **self.kwargs)
        self._sc.append(sc)
        self.sc = sc
        # cb = self.ax.figure.colorbar(sc, ax=self.ax)
        # cb.set_label(self.I)
        super(LiveScatter, self).start(doc)

    def update(self, x, y, I):

        # if one is None all are
        if self._minx is None:
            self._minx = self.xlim[0]
            self._maxx = self.xlim[1]
            self._miny = self.ylim[0]
            self._maxy = self.ylim[1]

        self._xdata.append(x)
        self._ydata.append(y)
        self._Idata.append(I)
        offsets = np.vstack([self._xdata, self._ydata]).T
        self.sc.set_offsets(offsets)
        self.sc.set_array(np.asarray(self._Idata))

        # if self.xlim is None:
 
        minx, maxx = np.minimum(x, self._minx), np.maximum(x, self._maxx)
        self.ax.set_xlim(minx, maxx)

        self._minx = minx
        self._maxx = maxx

        # if self.ylim is None:
        miny, maxy = np.minimum(y, self._miny), np.maximum(y, self._maxy)
        self.ax.set_ylim(miny, maxy)

        self._miny = miny
        self._maxy = maxy     

        if self.clim is None:
            clim = np.nanmin(self._Idata), np.nanmax(self._Idata)
            self.sc.set_clim(*clim)


class BestEffortCallbackHxn(BestEffortCallback):

       def descriptor(self, doc):
           
        self._descriptors[doc['uid']] = doc
        stream_name = doc.get('name', 'primary')  # fall back for old docs

        if stream_name not in self._stream_names_seen:
            self._stream_names_seen.add(stream_name)
            if self._table_enabled:
                print("New stream: {!r}".format(stream_name))

        columns = hinted_fields(doc)

        # ## This deals with old documents. ## #

        if stream_name == 'primary' and self._cleanup_motor_heuristic:
            # We stashed object names in self.dim_fields, which we now need to
            # look up the actual fields for.
            self._cleanup_motor_heuristic = False
            fixed_dim_fields = []
            for obj_name in self.dim_fields:
                # Special case: 'time' can be a dim_field, but it's not an
                # object name. Just add it directly to the list of fields.
                if obj_name == 'time':
                    fixed_dim_fields.append('time')
                    continue
                try:
                    fields = doc.get('hints', {}).get(obj_name, {})['fields']
                except KeyError:
                    fields = doc['object_keys'][obj_name]
                fixed_dim_fields.extend(fields)
            self.dim_fields = fixed_dim_fields

        # ## TABLE ## #

        if stream_name == self.dim_stream:
            # Ensure that no independent variables ('dimensions') are
            # duplicated here.
            columns = [c for c in columns if c not in self.all_dim_fields]

            if self._table_enabled:
                # plot everything, independent or dependent variables
                self._table = LiveTable(list(self.all_dim_fields) + columns)
                self._table('start', self._start_doc)
                self._table('descriptor', doc)

        # ## DECIDE WHICH KIND OF PLOT CAN BE USED ## #

        if not self._plots_enabled:
            return
        if stream_name in self.noplot_streams:
            return
        if not columns:
            return
        if ((self._start_doc.get('num_points') == 1) and
                (stream_name == self.dim_stream) and
                self.omit_single_point_plot):
            return

        # This is a heuristic approach until we think of how to hint this in a
        # generalizable way.
        if stream_name == self.dim_stream:
            dim_fields = self.dim_fields
        else:
            dim_fields = ['time']  # 'time' once LivePlot can do that

        # Create a figure or reuse an existing one.

        fig_name = '{} vs {}'.format(' '.join(sorted(columns)),
                                     ' '.join(sorted(dim_fields)))
        if self.overplot and len(dim_fields) == 1:
            # If any open figure matches 'figname {number}', use it. If there
            # are multiple, the most recently touched one will be used.
            pat1 = re.compile('^' + fig_name + '$')
            pat2 = re.compile('^' + fig_name + r' \d+$')
            for label in plt.get_figlabels():
                if pat1.match(label) or pat2.match(label):
                    fig_name = label
                    break
        else:
            if plt.fignum_exists(fig_name):
                # Generate a unique name by appending a number.
                for number in itertools.count(2):
                    new_name = '{} {}'.format(fig_name, number)
                    if not plt.fignum_exists(new_name):
                        fig_name = new_name
                        break
        ndims = len(dim_fields)
        if not 0 < ndims < 3:
            # we need 1 or 2 dims to do anything, do not make empty figures
            return

        if self._fig_factory:
            fig = self._fig_factory(fig_name)
        else:
            fig = plt.figure(fig_name)
        if not fig.axes:
            # This is apparently a fresh figure. Make axes.
            # The complexity here is due to making a shared x axis. This can be
            # simplified when Figure supports the `subplots` method in a future
            # release of matplotlib.
            fig.set_size_inches(6.4, min(950, len(columns) * 400) / fig.dpi)
            for i in range(len(columns)):
                if i == 0:
                    ax = fig.add_subplot(len(columns), 1, 1 + i)
                    if ndims == 1:
                        share_kwargs = {'sharex': ax}
                    elif ndims == 2:
                        share_kwargs = {'sharex': ax, 'sharey': ax}
                    else:
                        raise NotImplementedError("we now support 3D?!")
                else:
                    ax = fig.add_subplot(len(columns), 1, 1 + i,
                                         **share_kwargs)
        axes = fig.axes

        # ## LIVE PLOT AND PEAK ANALYSIS ## #

        if ndims == 1:
            self._live_plots[doc['uid']] = {}
            self._peak_stats[doc['uid']] = {}
            x_key, = dim_fields
            for y_key, ax in zip(columns, axes):
                dtype = doc['data_keys'][y_key]['dtype']
                if dtype not in ('number', 'integer'):
                    warn("Omitting {} from plot because dtype is {}"
                         "".format(y_key, dtype))
                    continue
                # Create an instance of LivePlot and an instance of PeakStats.
                live_plot = LivePlotPlusPeaks(y=y_key, x=x_key, ax=ax,
                                              peak_results=self.peaks)
                live_plot('start', self._start_doc)
                live_plot('descriptor', doc)
                peak_stats = PeakStats(x=x_key, y=y_key)
                peak_stats('start', self._start_doc)
                peak_stats('descriptor', doc)

                # Stash them in state.
                self._live_plots[doc['uid']][y_key] = live_plot
                self._peak_stats[doc['uid']][y_key] = peak_stats

            for ax in axes[:-1]:
                ax.set_xlabel('')
        elif ndims == 2:
            # Decide whether to use LiveGrid or LiveScatter. LiveScatter is the
            # safer one to use, so it is the fallback..
            gridding = self._start_doc.get('hints', {}).get('gridding')
            if gridding == 'rectilinear':
                self._live_grids[doc['uid']] = {}
                slow, fast = dim_fields
                try:
                    extents = self._start_doc['extents']
                    shape = self._start_doc['shape']
                except KeyError:
                    warn("Need both 'shape' and 'extents' in plan metadata to "
                         "create LiveGrid.")
                else:
                    data_range = np.array([float(np.diff(e)) for e in extents])
                    y_step, x_step = data_range / [max(1, s - 1) for s in shape]
                    adjusted_extent = [extents[1][0] - x_step / 2,
                                       extents[1][1] + x_step / 2,
                                       extents[0][0] - y_step / 2,
                                       extents[0][1] + y_step / 2]
                    for I_key, ax in zip(columns, axes):
                        # MAGIC NUMBERS based on what tacaswell thinks looks OK
                        data_aspect_ratio = np.abs(data_range[1]/data_range[0])
                        MAR = 2
                        if (1/MAR < data_aspect_ratio < MAR):
                            aspect = 'equal'
                            ax.set_aspect(aspect, adjustable='box-forced')
                        else:
                            aspect = 'auto'
                            ax.set_aspect(aspect, adjustable='datalim')

                        live_grid = LiveGrid(shape, I_key,
                                             xlabel=fast, ylabel=slow,
                                             extent=adjusted_extent,
                                             aspect=aspect,
                                             ax=ax)

                        live_grid('start', self._start_doc)
                        live_grid('descriptor', doc)
                        self._live_grids[doc['uid']][I_key] = live_grid
            else:
                self._live_scatters[doc['uid']] = {}
                x_key, y_key = dim_fields
                for I_key, ax in zip(columns, axes):
                    try:
                        extents = self._start_doc['extents']
                    except KeyError:
                        xlim = ylim = None
                    else:
                        xlim, ylim = extents
                    live_scatter = LiveScatterHxn(x_key, y_key, I_key,
                                               xlim=xlim, ylim=ylim,
                                               # Let clim autoscale.
                                               ax=ax)
                     
                    live_scatter('start', self._start_doc)
                    live_scatter('descriptor', doc)
                    self._live_scatters[doc['uid']][I_key] = live_scatter
        else:
            raise NotImplementedError("we do not support 3D+ in BEC yet "
                                      "(and it should have bailed above)")
        try:
            fig.tight_layout()
        except ValueError:
            pass


bec_hxn = BestEffortCallbackHxn()
# fermat_tst = bpp.subs_decorator(bec_fermat)(hxntools.scans.relative_fermat)

from IPython import get_ipython
user_ns = get_ipython().user_ns
ns = {}
ns['bec_hxn'] =  bec_hxn
user_ns.update(ns)

