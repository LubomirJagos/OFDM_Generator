#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: OFDM Tx
# Description: Example of an OFDM Transmitter
# Generated: Sat May 13 02:21:38 2017
##################################################

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print "Warning: failed to XInitThreads()"

from PyQt4 import Qt
from gnuradio import blocks
from gnuradio import digital
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import qtgui
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import numpy
import sip
import sys


class tx_ofdm(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "OFDM Tx")
        Qt.QWidget.__init__(self)
        self.setWindowTitle("OFDM Tx")
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except:
            pass
        self.top_scroll_layout = Qt.QVBoxLayout()
        self.setLayout(self.top_scroll_layout)
        self.top_scroll = Qt.QScrollArea()
        self.top_scroll.setFrameStyle(Qt.QFrame.NoFrame)
        self.top_scroll_layout.addWidget(self.top_scroll)
        self.top_scroll.setWidgetResizable(True)
        self.top_widget = Qt.QWidget()
        self.top_scroll.setWidget(self.top_widget)
        self.top_layout = Qt.QVBoxLayout(self.top_widget)
        self.top_grid_layout = Qt.QGridLayout()
        self.top_layout.addLayout(self.top_grid_layout)

        self.settings = Qt.QSettings("GNU Radio", "tx_ofdm")
        self.restoreGeometry(self.settings.value("geometry").toByteArray())

        ##################################################
        # Variables
        ##################################################
        self.occupied_carriers = occupied_carriers = (range(-102,-75) + range(-74,-39) + range(-38, -11) + range(-10, 0) + range(1, 11) + range(12, 39) + range(40,75) + range(76,103),)
        self.length_tag_key = length_tag_key = "packet_len"
        self.sync_word2 = sync_word2 = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,0,-1,0,-1,0,1,0,1,0,1,0,-1,0,-1,0,-1,0,1,0,-1,0,1,0,1,0,-1,0,-1,0,-1,0,-1,0,-1,0,1,0,-1,0,-1,0,1,0,-1,0,-1,0,1,0,-1,0,1,0,1,0,1,0,-1,0,1,0,-1,0,0,0,1,0,1,0,1,0,-1,0,1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,1,0,-1,0,1,0,1,0,1,0,-1,0,1,0,1,0,1,0,-1,0,1,0,1,0,1,0,1,0,1,0,-1,0,-1,0,1,0,-1,0,-1,0,-1,0,-1,0,1,0,1,0,0,0,-1,0,-1,0,1,0,-1,0,1,0,1,0,-1,0,-1,0,-1,0,1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,1,0,1,0,-1,0,1,0,-1,0,-1,0,1,0,1,0,-1,0,1,0,-1,0,1,0,1,0,1,0,1,0,1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
        self.sync_word1 = sync_word1 = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,1.4142,0,1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,1.4142,0,1.4142,0,1.4142,0,1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,1.4142,0,-1.4142,0,1.4142,0,1.4142,0,1.4142,0,1.4142,0,0,0,1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,1.4142,0,1.4142,0,1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,0,0,-1.4142,0,-1.4142,0,1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,1.4142,0,1.4142,0,-1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,1.4142,0,1.4142,0,1.4142,0,-1.4142,0,-1.4142,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
        self.samp_rate = samp_rate = 100000
        self.rolloff = rolloff = 0
        self.pilot_symbols = pilot_symbols = ((1,1,1,1,1,1,1,-1),)
        self.pilot_carriers = pilot_carriers = ((-103,-75,-39,-11,11,39,75,103),)
        self.payload_mod = payload_mod = digital.constellation_qpsk()
        self.packet_len = packet_len = 420
        self.header_mod = header_mod = digital.constellation_bpsk()
        self.hdr_format = hdr_format = digital.header_format_ofdm(occupied_carriers, 1, length_tag_key,)
        self.fft_len = fft_len = 256

        ##################################################
        # Blocks
        ##################################################
        self.qtgui_time_sink_x_0 = qtgui.time_sink_c(
        	1024, #size
        	samp_rate, #samp_rate
        	"", #name
        	1 #number of inputs
        )
        self.qtgui_time_sink_x_0.set_update_time(0.10)
        self.qtgui_time_sink_x_0.set_y_axis(-1, 1)
        
        self.qtgui_time_sink_x_0.set_y_label('Amplitude', "")
        
        self.qtgui_time_sink_x_0.enable_tags(-1, True)
        self.qtgui_time_sink_x_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, 0, "")
        self.qtgui_time_sink_x_0.enable_autoscale(True)
        self.qtgui_time_sink_x_0.enable_grid(False)
        self.qtgui_time_sink_x_0.enable_axis_labels(True)
        self.qtgui_time_sink_x_0.enable_control_panel(False)
        
        if not True:
          self.qtgui_time_sink_x_0.disable_legend()
        
        labels = ['', '', '', '', '',
                  '', '', '', '', '']
        widths = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        colors = ["blue", "red", "green", "black", "cyan",
                  "magenta", "yellow", "dark red", "dark green", "blue"]
        styles = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        markers = [-1, -1, -1, -1, -1,
                   -1, -1, -1, -1, -1]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0]
        
        for i in xrange(2*1):
            if len(labels[i]) == 0:
                if(i % 2 == 0):
                    self.qtgui_time_sink_x_0.set_line_label(i, "Re{{Data {0}}}".format(i/2))
                else:
                    self.qtgui_time_sink_x_0.set_line_label(i, "Im{{Data {0}}}".format(i/2))
            else:
                self.qtgui_time_sink_x_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_0.set_line_alpha(i, alphas[i])
        
        self._qtgui_time_sink_x_0_win = sip.wrapinstance(self.qtgui_time_sink_x_0.pyqwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_time_sink_x_0_win)
        self.digital_ofdm_tx_0 = digital.ofdm_tx(
        	  fft_len=fft_len, cp_len=fft_len/4,
        	  packet_length_tag_key='packet_len',
        	  occupied_carriers=occupied_carriers,
        	  pilot_carriers=pilot_carriers,
        	  pilot_symbols=pilot_symbols,
        	  sync_word1=sync_word1,
        	  sync_word2=sync_word2,
        	  bps_header=1,
        	  bps_payload=2,
        	  rolloff=0,
        	  debug_log=False,
        	  scramble_bits=False
        	 )
        self.digital_ofdm_rx_0 = digital.ofdm_rx(
        	  fft_len=fft_len, cp_len=fft_len/4,
        	  frame_length_tag_key='frame_'+'rx_len',
        	  packet_length_tag_key='rx_len',
        	  occupied_carriers=occupied_carriers,
        	  pilot_carriers=pilot_carriers,
        	  pilot_symbols=pilot_symbols,
        	  sync_word1=sync_word1,
        	  sync_word2=sync_word2,
        	  bps_header=1,
        	  bps_payload=2,
        	  debug_log=False,
        	  scramble_bits=False
        	 )
        self.blocks_throttle_0 = blocks.throttle(gr.sizeof_char*1, samp_rate,True)
        self.blocks_tag_debug_0 = blocks.tag_debug(gr.sizeof_char*1, '', ""); self.blocks_tag_debug_0.set_display(True)
        self.blocks_stream_to_tagged_stream_0 = blocks.stream_to_tagged_stream(gr.sizeof_char, 1, packet_len, "packet_len")
        self.blocks_file_sink_0_0 = blocks.file_sink(gr.sizeof_gr_complex*1, 'C:\\Users\\Miska\\Lubo_Diplomka\\_Matlab_ostra_verzia\\gnuradio_ofdm_headerProblem\\ofdm_sig_out.txt', False)
        self.blocks_file_sink_0_0.set_unbuffered(False)
        self.blocks_file_sink_0 = blocks.file_sink(gr.sizeof_char*1, 'C:\\Users\\Miska\\Lubo_Diplomka\\_Matlab_ostra_verzia\\gnuradio_ofdm_headerProblem\\ofdm_sig_in.txt', False)
        self.blocks_file_sink_0.set_unbuffered(False)
        self.analog_random_source_x_0 = blocks.vector_source_b(map(int, numpy.random.randint(0, 255, 1000)), True)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.analog_random_source_x_0, 0), (self.blocks_stream_to_tagged_stream_0, 0))    
        self.connect((self.blocks_stream_to_tagged_stream_0, 0), (self.blocks_throttle_0, 0))    
        self.connect((self.blocks_throttle_0, 0), (self.blocks_file_sink_0, 0))    
        self.connect((self.blocks_throttle_0, 0), (self.digital_ofdm_tx_0, 0))    
        self.connect((self.digital_ofdm_rx_0, 0), (self.blocks_tag_debug_0, 0))    
        self.connect((self.digital_ofdm_tx_0, 0), (self.blocks_file_sink_0_0, 0))    
        self.connect((self.digital_ofdm_tx_0, 0), (self.digital_ofdm_rx_0, 0))    
        self.connect((self.digital_ofdm_tx_0, 0), (self.qtgui_time_sink_x_0, 0))    

    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "tx_ofdm")
        self.settings.setValue("geometry", self.saveGeometry())
        event.accept()

    def get_occupied_carriers(self):
        return self.occupied_carriers

    def set_occupied_carriers(self, occupied_carriers):
        self.occupied_carriers = occupied_carriers
        self.set_hdr_format(digital.header_format_ofdm(self.occupied_carriers, 1, self.length_tag_key,))

    def get_length_tag_key(self):
        return self.length_tag_key

    def set_length_tag_key(self, length_tag_key):
        self.length_tag_key = length_tag_key
        self.set_hdr_format(digital.header_format_ofdm(self.occupied_carriers, 1, self.length_tag_key,))

    def get_sync_word2(self):
        return self.sync_word2

    def set_sync_word2(self, sync_word2):
        self.sync_word2 = sync_word2

    def get_sync_word1(self):
        return self.sync_word1

    def set_sync_word1(self, sync_word1):
        self.sync_word1 = sync_word1

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.qtgui_time_sink_x_0.set_samp_rate(self.samp_rate)
        self.blocks_throttle_0.set_sample_rate(self.samp_rate)

    def get_rolloff(self):
        return self.rolloff

    def set_rolloff(self, rolloff):
        self.rolloff = rolloff

    def get_pilot_symbols(self):
        return self.pilot_symbols

    def set_pilot_symbols(self, pilot_symbols):
        self.pilot_symbols = pilot_symbols

    def get_pilot_carriers(self):
        return self.pilot_carriers

    def set_pilot_carriers(self, pilot_carriers):
        self.pilot_carriers = pilot_carriers

    def get_payload_mod(self):
        return self.payload_mod

    def set_payload_mod(self, payload_mod):
        self.payload_mod = payload_mod

    def get_packet_len(self):
        return self.packet_len

    def set_packet_len(self, packet_len):
        self.packet_len = packet_len
        self.blocks_stream_to_tagged_stream_0.set_packet_len(self.packet_len)
        self.blocks_stream_to_tagged_stream_0.set_packet_len_pmt(self.packet_len)

    def get_header_mod(self):
        return self.header_mod

    def set_header_mod(self, header_mod):
        self.header_mod = header_mod

    def get_hdr_format(self):
        return self.hdr_format

    def set_hdr_format(self, hdr_format):
        self.hdr_format = hdr_format

    def get_fft_len(self):
        return self.fft_len

    def set_fft_len(self, fft_len):
        self.fft_len = fft_len


def main(top_block_cls=tx_ofdm, options=None):

    from distutils.version import StrictVersion
    if StrictVersion(Qt.qVersion()) >= StrictVersion("4.5.0"):
        style = gr.prefs().get_string('qtgui', 'style', 'raster')
        Qt.QApplication.setGraphicsSystem(style)
    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls()
    tb.start()
    tb.show()

    def quitting():
        tb.stop()
        tb.wait()
    qapp.connect(qapp, Qt.SIGNAL("aboutToQuit()"), quitting)
    qapp.exec_()


if __name__ == '__main__':
    main()
