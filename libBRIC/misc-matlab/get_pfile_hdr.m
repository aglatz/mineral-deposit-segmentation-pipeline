function [hdr] = get_pfile_hdr(gehdr)

hdr = struct;

hdr.data_offset = gehdr.rdb.off_data;
hdr.point_size = gehdr.rdb.point_size;

hdr.nsamp = gehdr.rdb.da_xres;
hdr.nframe = gehdr.rdb.da_yres;

start_recv = gehdr.rdb.dab_start_rcv(1);
stop_recv = gehdr.rdb.dab_stop_rcv(1);
hdr.nrecv = (stop_recv - start_recv) + 1;
 
hdr.necho = gehdr.rdb.nechoes;
hdr.nslice = gehdr.rdb.nslices/gehdr.rdb.npasses;