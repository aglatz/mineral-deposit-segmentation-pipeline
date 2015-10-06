function hdr = read_gehdr(id)
%READ_GEHDR  Read the header of a GE raw data file into a structure
%
%  usage: hdr = read_gehdr(id)
%
%  id:    integer file identifier obtained from fopen()
%
%  hdr:   a structure containing the following elements:
%
%         hdr.rdb - a structure equivilant to GE's rdb_hdr_rec
%         hdr.image - a structure equivilant to GE's rdb_hdr_image
%         hdr.series - a structure equivilant to GE's rdb_hdr_series
%         hdr.exam - a structure equivilant to GE's rdb_hdr_exam
%         hdr.data_acq_tab - a structure array equivilant to GE's rdb_hdr_data_acq_tab
%
%         and these elements that are straight storage dumps of their respective
%         GE structures:
%
%         hdr.per_pass
%         hdr.unlock_raw
%         hdr.nex_tab
%         hdr.nex_abort_tab
%         hdr.tool
%
%  e.g.:  >> id = fopen('P12345.7', 'r', 'l');
%         >> hdr = read_gehdr(id);
%
%  note:  Raw file must have RDB_RDBM_REVISION = 15.001F

%  This script was generated from GE source using 'gehdr2matlab' written by DB Clayton


fseek(id, 0, -1);  % go to start of raw data file
hdr.rdb.rdbm_rev = fread(id, 1, 'float');  
hdr.rdb.run_int = fread(id, 1, 'long');  % Rdy pkt Run Number 
hdr.rdb.scan_seq = fread(id, 1, 'short');  % Rdy pkt Sequence Number 
hdr.rdb.run_char = fread(id, 6, '*char')';  % Rdy pkt Run no in char 
hdr.rdb.scan_date = fread(id, 10, '*char')';  % 
hdr.rdb.scan_time = fread(id, 8, '*char')';  % 
hdr.rdb.logo = fread(id, 10, '*char')';  % rdbm used to verify file 
hdr.rdb.file_contents = fread(id, 1, 'short');  % Data type 0=emp 1=nrec 2=rw 0, 1, 2 
hdr.rdb.lock_mode = fread(id, 1, 'short');  % unused 
hdr.rdb.dacq_ctrl = fread(id, 1, 'short');  % rhdacqctrl bit mask 15 bits 
hdr.rdb.recon_ctrl = fread(id, 1, 'short');  % rhrcctrl bit mask 15 bits 
hdr.rdb.exec_ctrl = fread(id, 1, 'short');  % rhexecctrl bit mask 15 bits 
hdr.rdb.scan_type = fread(id, 1, 'short');  % bit mask 15 bits 
hdr.rdb.data_collect_type = fread(id, 1, 'short');  % rhtype bit mask 15 bits 
hdr.rdb.data_format = fread(id, 1, 'short');  % rhformat bit mask 15 bits 
hdr.rdb.recon = fread(id, 1, 'short');  % rhrecon proc-a-son recon 0 - 100 
hdr.rdb.datacq = fread(id, 1, 'short');  % rhdatacq proc-a-son dacq 
hdr.rdb.npasses = fread(id, 1, 'short');  % rhnpasses passes for a scan 0 - 256 
hdr.rdb.npomp = fread(id, 1, 'short');  % rhnpomp pomp group slices 1,2 
hdr.rdb.nslices = fread(id, 1, 'ushort');  % rhnslices slices in a pass 0 - 256 
hdr.rdb.nechoes = fread(id, 1, 'short');  % rhnecho echoes of a slice 1 - 32 
hdr.rdb.navs = fread(id, 1, 'short');  % rhnavs num of excitiations 1 - 32727 
hdr.rdb.nframes = fread(id, 1, 'short');  % rhnframes yres 0 - 1024 
hdr.rdb.baseline_views = fread(id, 1, 'short');  % rhbline baselines 0 - 1028 
hdr.rdb.hnover = fread(id, 1, 'short');  % rhhnover overscans 0 - 1024 
hdr.rdb.frame_size = fread(id, 1, 'ushort');  % rhfrsize xres 0 - 32768 
hdr.rdb.point_size = fread(id, 1, 'short');  % rhptsize 2 - 4 
hdr.rdb.vquant = fread(id, 1, 'short');  % rhvquant 3d volumes 1 
hdr.rdb.cheart = fread(id, 1, 'short');  % RX Cine heart phases 1 - 32 
hdr.rdb.ctr = fread(id, 1, 'float');  % RX Cine TR in sec 0 - 3.40282e38
hdr.rdb.ctrr = fread(id, 1, 'float');  % RX Cine RR in sec 0 - 30.0 
hdr.rdb.initpass = fread(id, 1, 'short');  % rhinitpass allocate passes 0 - 32767 
hdr.rdb.incrpass = fread(id, 1, 'short');  % rhincrpass tps autopauses 0 - 32767 
hdr.rdb.method_ctrl = fread(id, 1, 'short');  % rhmethod 0=recon, 1=psd 0, 1 
hdr.rdb.da_xres = fread(id, 1, 'ushort');  % rhdaxres 0 - 32768 
hdr.rdb.da_yres = fread(id, 1, 'short');  % rhdayres 0 - 2049 
hdr.rdb.rc_xres = fread(id, 1, 'short');  % rhrcxres 0 - 1024 
hdr.rdb.rc_yres = fread(id, 1, 'short');  % rhrcyres 0 - 1024 
hdr.rdb.im_size = fread(id, 1, 'short');  % rhimsize 0 - 512 
hdr.rdb.rc_zres = fread(id, 1, 'long');  % power of 2 > rhnslices 0 - 128 
hdr.rdb.raw_pass_size = fread(id, 1, 'ulong');  % rhrawsize 0 - 2147483647
hdr.rdb.sspsave = fread(id, 1, 'ulong');  % rhsspsave 0 - 2147483647
hdr.rdb.udasave = fread(id, 1, 'ulong');  % rhudasave 0 - 2147483647
hdr.rdb.fermi_radius = fread(id, 1, 'float');  % rhfermr fermi radius 0 - 3.40282e38
hdr.rdb.fermi_width = fread(id, 1, 'float');  % rhfermw fermi width 0 - 3.40282e38
hdr.rdb.fermi_ecc = fread(id, 1, 'float');  % rhferme fermi excentiricty 0 - 3.40282e38
hdr.rdb.clip_min = fread(id, 1, 'float');  % rhclipmin 4x IP limit +-16383 
hdr.rdb.clip_max = fread(id, 1, 'float');  % rhclipmax 4x IP limit +-16383 
hdr.rdb.default_offset = fread(id, 1, 'float');  % rhdoffset default offset = 0 +-3.40282e38 
hdr.rdb.xoff = fread(id, 1, 'float');  % rhxoff scroll img in x +-256 
hdr.rdb.yoff = fread(id, 1, 'float');  % rhyoff scroll img in y +-256 
hdr.rdb.nwin = fread(id, 1, 'float');  % rhnwin hecho window width 0 - 256 
hdr.rdb.ntran = fread(id, 1, 'float');  % rhntran hecho trans width 0 - 256 
hdr.rdb.scalei = fread(id, 1, 'float');  % PS rhscalei +-3.40282e38 
hdr.rdb.scaleq = fread(id, 1, 'float');  % PS rhscaleq def = 0 +-3.40282e38 
hdr.rdb.rotation = fread(id, 1, 'short');  % RX 0 90 180 270 deg 0 - 3 
hdr.rdb.transpose = fread(id, 1, 'short');  % RX 0, 1 n / y transpose 0 - 1
hdr.rdb.kissoff_views = fread(id, 1, 'short');  % rhblank zero image views 0 - 512 
hdr.rdb.slblank = fread(id, 1, 'short');  % rhslblank slice blank 3d 0 - 128 
hdr.rdb.gradcoil = fread(id, 1, 'short');  % RX 0=off 1=Schnk 2=Rmr 0 - 2 
hdr.rdb.ddaover = fread(id, 1, 'short');  % rhddaover unused 
hdr.rdb.sarr = fread(id, 1, 'short');  % SARR bit mask 15 bits 
hdr.rdb.fd_tr = fread(id, 1, 'short');  % SARR feeder timing info 
hdr.rdb.fd_te = fread(id, 1, 'short');  % SARR feeder timing info 
hdr.rdb.fd_ctrl = fread(id, 1, 'short');  % SARR control of feeder 
hdr.rdb.algor_num = fread(id, 1, 'short');  % SARR df decimation ratio 
hdr.rdb.fd_df_dec = fread(id, 1, 'short');  % SARR which feeder algor 
buff = fread(id, 8, 'short');  % kluge for type RDB_MULTI_RCV_TYPE
hdr.rdb.dab_start_rcv = buff(1:2:end);  % kluge for type RDB_MULTI_RCV_TYPE
hdr.rdb.dab_stop_rcv = buff(2:2:end);  % kluge for type RDB_MULTI_RCV_TYPE
hdr.rdb.user0 = fread(id, 1, 'float');  % rhuser0 +-3.40282e38 
hdr.rdb.user1 = fread(id, 1, 'float');  % rhuser1 +-3.40282e38 
hdr.rdb.user2 = fread(id, 1, 'float');  % rhuser2 +-3.40282e38 
hdr.rdb.user3 = fread(id, 1, 'float');  % rhuser3 +-3.40282e38 
hdr.rdb.user4 = fread(id, 1, 'float');  % rhuser4 +-3.40282e38 
hdr.rdb.user5 = fread(id, 1, 'float');  % rhuser5 +-3.40282e38 
hdr.rdb.user6 = fread(id, 1, 'float');  % rhuser6 +-3.40282e38 
hdr.rdb.user7 = fread(id, 1, 'float');  % rhuser7 +-3.40282e38 
hdr.rdb.user8 = fread(id, 1, 'float');  % rhuser8 +-3.40282e38 
hdr.rdb.user9 = fread(id, 1, 'float');  % rhuser9 +-3.40282e38 
hdr.rdb.user10 = fread(id, 1, 'float');  % rhuser10 +-3.40282e38 
hdr.rdb.user11 = fread(id, 1, 'float');  % rhuser11 +-3.40282e38 
hdr.rdb.user12 = fread(id, 1, 'float');  % rhuser12 +-3.40282e38 
hdr.rdb.user13 = fread(id, 1, 'float');  % rhuser13 +-3.40282e38 
hdr.rdb.user14 = fread(id, 1, 'float');  % rhuser14 +-3.40282e38 
hdr.rdb.user15 = fread(id, 1, 'float');  % rhuser15 +-3.40282e38 
hdr.rdb.user16 = fread(id, 1, 'float');  % rhuser16 +-3.40282e38 
hdr.rdb.user17 = fread(id, 1, 'float');  % rhuser17 +-3.40282e38 
hdr.rdb.user18 = fread(id, 1, 'float');  % rhuser18 +-3.40282e38 
hdr.rdb.user19 = fread(id, 1, 'float');  % rhuser19 +-3.40282e38 
hdr.rdb.v_type = fread(id, 1, 'long');  % rhvtype bit mask 31 bits 
hdr.rdb.v_coefxa = fread(id, 1, 'float');  % RX x flow direction control 0 - 4 
hdr.rdb.v_coefxb = fread(id, 1, 'float');  % RX x flow direction control 0 - 4 
hdr.rdb.v_coefxc = fread(id, 1, 'float');  % RX x flow direction control 0 - 4 
hdr.rdb.v_coefxd = fread(id, 1, 'float');  % RX x flow direction control 0 - 4 
hdr.rdb.v_coefya = fread(id, 1, 'float');  % RX y flow direction control 0 - 4 
hdr.rdb.v_coefyb = fread(id, 1, 'float');  % RX y flow direction control 0 - 4 
hdr.rdb.v_coefyc = fread(id, 1, 'float');  % RX y flow direction control 0 - 4 
hdr.rdb.v_coefyd = fread(id, 1, 'float');  % RX y flow direction control 0 - 4 
hdr.rdb.v_coefza = fread(id, 1, 'float');  % RX z flow direction control 0 - 4 
hdr.rdb.v_coefzb = fread(id, 1, 'float');  % RX z flow direction control 0 - 4 
hdr.rdb.v_coefzc = fread(id, 1, 'float');  % RX z flow direction control 0 - 4 
hdr.rdb.v_coefzd = fread(id, 1, 'float');  % RX z flow direction control 0 - 4 
hdr.rdb.vm_coef1 = fread(id, 1, 'float');  % RX weight for mag image 1 0 - 1 
hdr.rdb.vm_coef2 = fread(id, 1, 'float');  % RX weight for mag image 2 0 - 1 
hdr.rdb.vm_coef3 = fread(id, 1, 'float');  % RX weight for mag image 3 0 - 1 
hdr.rdb.vm_coef4 = fread(id, 1, 'float');  % RX weight for mag image 4 0 - 1 
hdr.rdb.v_venc = fread(id, 1, 'float');  % RX vel encodeing cm / sec 0.001 - 5000 
hdr.rdb.spectral_width = fread(id, 1, 'float');  % specwidth filter width kHz 500 - 3355432 
hdr.rdb.csi_dims = fread(id, 1, 'short');  % spectro 
hdr.rdb.xcsi = fread(id, 1, 'short');  % rhspecrescsix 2 - 64 
hdr.rdb.ycsi = fread(id, 1, 'short');  % rhspecrescsiy 2 - 64 
hdr.rdb.zcsi = fread(id, 1, 'short');  % spectro 
hdr.rdb.roilenx = fread(id, 1, 'float');  % RX x csi volume dimension 
hdr.rdb.roileny = fread(id, 1, 'float');  % RX y csi volume dimension 
hdr.rdb.roilenz = fread(id, 1, 'float');  % RX z csi volume dimension 
hdr.rdb.roilocx = fread(id, 1, 'float');  % RX x csi volume center 
hdr.rdb.roilocy = fread(id, 1, 'float');  % RX y csi volume center 
hdr.rdb.roilocz = fread(id, 1, 'float');  % RX z csi volume center 
hdr.rdb.numdwell = fread(id, 1, 'float');  % specdwells 0 - 3.40282e38
hdr.rdb.ps_command = fread(id, 1, 'long');  % PS internal use only 
hdr.rdb.ps_mps_r1 = fread(id, 1, 'long');  % PS MPS R1 setting 1 - 7 
hdr.rdb.ps_mps_r2 = fread(id, 1, 'long');  % PS MPS R2 setting 1 - 30 
hdr.rdb.ps_mps_tg = fread(id, 1, 'long');  % PS MPS Transmit gain setting 0 - 200
hdr.rdb.ps_mps_freq = fread(id, 1, 'long');  % PS MPS Center frequency hz +-3.40282e38 
hdr.rdb.ps_aps_r1 = fread(id, 1, 'long');  % PS APS R1 setting 1 - 7 
hdr.rdb.ps_aps_r2 = fread(id, 1, 'long');  % PS APS R2 setting 1 - 30 
hdr.rdb.ps_aps_tg = fread(id, 1, 'long');  % PS APS Transmit gain setting 0 - 200
hdr.rdb.ps_aps_freq = fread(id, 1, 'long');  % PS APS Center frequency hz +-3.40282e38 
hdr.rdb.ps_scalei = fread(id, 1, 'float');  % PS rational scaling +-3.40282e38 
hdr.rdb.ps_scaleq = fread(id, 1, 'float');  % PS unused 
hdr.rdb.ps_snr_warning = fread(id, 1, 'long');  % PS noise test 0=16 1=32 bits 0, 1 
hdr.rdb.ps_aps_or_mps = fread(id, 1, 'long');  % PS prescan order logic 0 - 5 
hdr.rdb.ps_mps_bitmap = fread(id, 1, 'long');  % PS bit mask 4 bits
hdr.rdb.ps_powerspec = fread(id, 256, '*char')';  % PS 
hdr.rdb.ps_filler1 = fread(id, 1, 'long');  % PS filler 
hdr.rdb.ps_filler2 = fread(id, 1, 'long');  % PS filler 
hdr.rdb.obsolete1 = fread(id, 16, 'float');  % PS mean noise each receiver +-3.40282e38 
hdr.rdb.obsolete2 = fread(id, 16, 'float');  % PS noise calc for muti rec +-3.40282e38 
hdr.rdb.halfecho = fread(id, 1, 'short');  % spectro full, half echo 0, 1 
hdr.rdb.im_size_y = fread(id, 1, 'short');  % rh???? 0 - 512 
hdr.rdb.data_collect_type1 = fread(id, 1, 'long');  % rh???? bit mask 31 bits 
hdr.rdb.freq_scale = fread(id, 1, 'float');  % rh???? freq k-space step +-3.40282e38 
hdr.rdb.phase_scale = fread(id, 1, 'float');  % rh???? freq k-space step +-3.40282e38 
hdr.rdb.ovl = fread(id, 1, 'short');  % rhovl - overlaps for MOTSA 
hdr.rdb.pclin = fread(id, 1, 'short');  % Linear Corr. 0:off, 1:linear, 2:polynomial 
hdr.rdb.pclinnpts = fread(id, 1, 'short');  % fit number of points 
hdr.rdb.pclinorder = fread(id, 1, 'short');  % fit order 
hdr.rdb.pclinavg = fread(id, 1, 'short');  % linear phase corr avg 0:off, 1:on 
hdr.rdb.pccon = fread(id, 1, 'short');  % Const Corr. 0:off, 1:Ky spec., 2:polyfit(2/ilv), 3:polyfit(1/ilv) 
hdr.rdb.pcconnpts = fread(id, 1, 'short');  % fit number of points 
hdr.rdb.pcconorder = fread(id, 1, 'short');  % fit order 
hdr.rdb.pcextcorr = fread(id, 1, 'short');  % external correction file 0:don't use, 1: use 
hdr.rdb.pcgraph = fread(id, 1, 'short');  % Phase Correction coef. image 0:off, 1:linear & constant 
hdr.rdb.pcileave = fread(id, 1, 'short');  % Interleaves to use for correction: 0=all, 1=only first 
hdr.rdb.hdbestky = fread(id, 1, 'short');  % bestky view for fractional Ky scan 
hdr.rdb.pcctrl = fread(id, 1, 'short');  % phase correction research control 
hdr.rdb.pcthrespts = fread(id, 1, 'short');  % 2..512 adjacent points 
hdr.rdb.pcdiscbeg = fread(id, 1, 'short');  % 0..512 beginning point to discard 
hdr.rdb.pcdiscmid = fread(id, 1, 'short');  % 0..512 middle point to discard 
hdr.rdb.pcdiscend = fread(id, 1, 'short');  % 0..512 ending point to discard 
hdr.rdb.pcthrespct = fread(id, 1, 'short');  % Threshold percentage 
hdr.rdb.pcspacial = fread(id, 1, 'short');  % Spacial best ref scan index 0..512 
hdr.rdb.pctemporal = fread(id, 1, 'short');  % Temporal best ref scan index 0..512 
hdr.rdb.pcspare = fread(id, 1, 'short');  % spare for phase correction 
hdr.rdb.ileaves = fread(id, 1, 'short');  % Number of interleaves 
hdr.rdb.kydir = fread(id, 1, 'short');  % Ky traversal dircetion 0: top-down, 1:center out 
hdr.rdb.alt = fread(id, 1, 'short');  % Alt read sign 0=no, 1=odd/even, 2=pairs 
hdr.rdb.reps = fread(id, 1, 'short');  % Number of scan repetitions 
hdr.rdb.ref = fread(id, 1, 'short');  % Ref Scan 0: off 1: on 
hdr.rdb.pcconnorm = fread(id, 1, 'float');  % Constant S term normalization factor 
hdr.rdb.pcconfitwt = fread(id, 1, 'float');  % Constant polyfit weighting factor 
hdr.rdb.pclinnorm = fread(id, 1, 'float');  % Linear S term normalization factor 
hdr.rdb.pclinfitwt = fread(id, 1, 'float');  % Linear polyfit weighting factor 
hdr.rdb.pcbestky = fread(id, 1, 'float');  % Best Ky location 
hdr.rdb.vrgf = fread(id, 1, 'long');  % control word for VRG filter 
hdr.rdb.vrgfxres = fread(id, 1, 'long');  % control word for VRGF final x resolution 
hdr.rdb.bp_corr = fread(id, 1, 'long');  % control word for bandpass asymmetry 
hdr.rdb.recv_freq_s = fread(id, 1, 'float');  % starting frequency (+62.5) 
hdr.rdb.recv_freq_e = fread(id, 1, 'float');  % ending frequency (-62.5) 
hdr.rdb.hniter = fread(id, 1, 'long');  % Selects the number of (continued...)
hdr.rdb.fast_rec = fread(id, 1, 'long');  % Added for homodyne II, tells if (continued...)
hdr.rdb.refframes = fread(id, 1, 'long');  % total # of frames for ref scan 
hdr.rdb.refframep = fread(id, 1, 'long');  % # of frames per pass for a ref scan 
hdr.rdb.scnframe = fread(id, 1, 'long');  % total # of frames for a entire scan 
hdr.rdb.pasframe = fread(id, 1, 'long');  % # of frames per pass 
hdr.rdb.user_usage_tag = fread(id, 1, 'ulong');  % for spectro 
hdr.rdb.user_fill_mapMSW = fread(id, 1, 'ulong');  % for spectro 
hdr.rdb.user_fill_mapLSW = fread(id, 1, 'ulong');  % for Spectro 
hdr.rdb.user20 = fread(id, 1, 'float');  % all following usercv are for spectro 
hdr.rdb.user21 = fread(id, 1, 'float');  
hdr.rdb.user22 = fread(id, 1, 'float');  
hdr.rdb.user23 = fread(id, 1, 'float');  
hdr.rdb.user24 = fread(id, 1, 'float');  
hdr.rdb.user25 = fread(id, 1, 'float');  
hdr.rdb.user26 = fread(id, 1, 'float');  
hdr.rdb.user27 = fread(id, 1, 'float');  
hdr.rdb.user28 = fread(id, 1, 'float');  
hdr.rdb.user29 = fread(id, 1, 'float');  
hdr.rdb.user30 = fread(id, 1, 'float');  
hdr.rdb.user31 = fread(id, 1, 'float');  
hdr.rdb.user32 = fread(id, 1, 'float');  
hdr.rdb.user33 = fread(id, 1, 'float');  
hdr.rdb.user34 = fread(id, 1, 'float');  
hdr.rdb.user35 = fread(id, 1, 'float');  
hdr.rdb.user36 = fread(id, 1, 'float');  
hdr.rdb.user37 = fread(id, 1, 'float');  
hdr.rdb.user38 = fread(id, 1, 'float');  
hdr.rdb.user39 = fread(id, 1, 'float');  
hdr.rdb.user40 = fread(id, 1, 'float');  
hdr.rdb.user41 = fread(id, 1, 'float');  
hdr.rdb.user42 = fread(id, 1, 'float');  
hdr.rdb.user43 = fread(id, 1, 'float');  
hdr.rdb.user44 = fread(id, 1, 'float');  
hdr.rdb.user45 = fread(id, 1, 'float');  
hdr.rdb.user46 = fread(id, 1, 'float');  
hdr.rdb.user47 = fread(id, 1, 'float');  
hdr.rdb.user48 = fread(id, 1, 'float');  
hdr.rdb.pcfitorig = fread(id, 1, 'short');  % Adjust view indexes if set so bestky view = 0 
hdr.rdb.pcshotfirst = fread(id, 1, 'short');  % First view within an echo group used for fit 
hdr.rdb.pcshotlast = fread(id, 1, 'short');  % Last view within an echo group used for fit 
hdr.rdb.pcmultegrp = fread(id, 1, 'short');  % If = 1, force pts from other egrps to be used 
hdr.rdb.pclinfix = fread(id, 1, 'short');  % If = 2, force slope to be set to pclinslope 
hdr.rdb.pcconfix = fread(id, 1, 'short');  % If = 2, force slope to be set to pcconslope 
hdr.rdb.pclinslope = fread(id, 1, 'float');  % Value to set lin slope to if forced 
hdr.rdb.pcconslope = fread(id, 1, 'float');  % Value to set con slope to if forced 
hdr.rdb.pccoil = fread(id, 1, 'short');  % If 1,2,3,4, use that coil's results for all 
hdr.rdb.vvsmode = fread(id, 1, 'short');  % Variable view sharing mode 
hdr.rdb.vvsaimgs = fread(id, 1, 'short');  % number of original images 
hdr.rdb.vvstr = fread(id, 1, 'short');  % TR in microseconds 
hdr.rdb.vvsgender = fread(id, 1, 'short');  % gender: male or female 
hdr.rdb.zip_factor = fread(id, 1, 'short');  % Slice ZIP factor: 0=OFF, 2, or 4 
hdr.rdb.maxcoef1a = fread(id, 1, 'float');  % Coefficient A for flow image 1 
hdr.rdb.maxcoef1b = fread(id, 1, 'float');  % Coefficient B for flow image 1 
hdr.rdb.maxcoef1c = fread(id, 1, 'float');  % Coefficient C for flow image 1 
hdr.rdb.maxcoef1d = fread(id, 1, 'float');  % Coefficient D for flow image 1 
hdr.rdb.maxcoef2a = fread(id, 1, 'float');  % Coefficient A for flow image 2 
hdr.rdb.maxcoef2b = fread(id, 1, 'float');  % Coefficient B for flow image 2 
hdr.rdb.maxcoef2c = fread(id, 1, 'float');  % Coefficient C for flow image 2 
hdr.rdb.maxcoef2d = fread(id, 1, 'float');  % Coefficient D for flow image 2 
hdr.rdb.maxcoef3a = fread(id, 1, 'float');  % Coefficient A for flow image 3 
hdr.rdb.maxcoef3b = fread(id, 1, 'float');  % Coefficient B for flow image 3 
hdr.rdb.maxcoef3c = fread(id, 1, 'float');  % Coefficient C for flow image 3 
hdr.rdb.maxcoef3d = fread(id, 1, 'float');  % Coefficient D for flow image 3 
hdr.rdb.ut_ctrl = fread(id, 1, 'long');  % System utility control variable 
hdr.rdb.dp_type = fread(id, 1, 'short');  % EPI II diffusion control cv 
hdr.rdb.arw = fread(id, 1, 'short');  % Arrhythmia rejection window(percentage:1-100)
hdr.rdb.vps = fread(id, 1, 'short');  % View Per Segment for FastCine 
hdr.rdb.mcReconEnable = fread(id, 1, 'short');  % N-Coil recon map 
hdr.rdb.fov = fread(id, 1, 'float');  % Auto-NCoil 
hdr.rdb.te = fread(id, 1, 'long');  % TE for first echo 
hdr.rdb.te2 = fread(id, 1, 'long');  % TE for second and later echoes 
hdr.rdb.dfmrbw = fread(id, 1, 'float');  % BW for navigator frames 
hdr.rdb.dfmctrl = fread(id, 1, 'long');  % Control flag for dfm (0=off, other=on)
hdr.rdb.raw_nex = fread(id, 1, 'long');  % Uncombined NEX at start of recon 
hdr.rdb.navs_per_pass = fread(id, 1, 'long');  % Max. navigator frames in a pass 
hdr.rdb.dfmxres = fread(id, 1, 'long');  % xres of navigator frames 
hdr.rdb.dfmptsize = fread(id, 1, 'long');  % point size of navigator frames 
hdr.rdb.navs_per_view = fread(id, 1, 'long');  % Num. navigators per frame (tag table) 
hdr.rdb.dfmdebug = fread(id, 1, 'long');  % control flag for dfm debug 
hdr.rdb.dfmthreshold = fread(id, 1, 'float');  % threshold for navigator correction 
hdr.rdb.grid_control = fread(id, 1, 'short');  % bit settings controlling gridding 
hdr.rdb.b0map = fread(id, 1, 'short');  % B0 map enable and map size 
hdr.rdb.grid_tediff = fread(id, 1, 'short');  % TE difference between b0 map arms 
hdr.rdb.grid_motion_comp = fread(id, 1, 'short');  % flag to apply motion compensation 
hdr.rdb.grid_radius_a = fread(id, 1, 'float');  % variable density transition 
hdr.rdb.grid_radius_b = fread(id, 1, 'float');  % variable density transition 
hdr.rdb.grid_max_gradient = fread(id, 1, 'float');  % Max gradient amplitude 
hdr.rdb.grid_max_slew = fread(id, 1, 'float');  % Max slew rate 
hdr.rdb.grid_scan_fov = fread(id, 1, 'float');  % Rx scan field of view 
hdr.rdb.grid_a2d_time = fread(id, 1, 'float');  % A to D sample time microsecs 
hdr.rdb.grid_density_factor = fread(id, 1, 'float');  % change factor for variable density 
hdr.rdb.grid_display_fov = fread(id, 1, 'float');  % Rx display field of view 
hdr.rdb.fatwater = fread(id, 1, 'short');  % for Fat and Water Dual Recon 
hdr.rdb.fiestamlf = fread(id, 1, 'short');  % MFO FIESTA recon control bit 16bits 
hdr.rdb.app = fread(id, 1, 'short');  % Auto Post-Processing opcode 
hdr.rdb.rhncoilsel = fread(id, 1, 'short');  % Auto-Ncoil 
hdr.rdb.rhncoillimit = fread(id, 1, 'short');  % Auto-Ncoil 
hdr.rdb.app_option = fread(id, 1, 'short');  % Auto Post_processing options 
hdr.rdb.grad_mode = fread(id, 1, 'short');  % Gradient mode in Gemini project 
hdr.rdb.pfile_passes = fread(id, 1, 'short');  % Num passes stored in a multi-pass Pfile (0 means 1 pass) 
hdr.rdb.asset = fread(id, 1, 'int');  
hdr.rdb.asset_calthresh = fread(id, 1, 'int');  
hdr.rdb.asset_R = fread(id, 1, 'float');  
hdr.rdb.coilno = fread(id, 1, 'int');  
hdr.rdb.asset_phases = fread(id, 1, 'int');  
hdr.rdb.scancent = fread(id, 1, 'float');  % Table position 
hdr.rdb.position = fread(id, 1, 'int');  % Patient position 
hdr.rdb.entry = fread(id, 1, 'int');  % Patient entry 
hdr.rdb.lmhor = fread(id, 1, 'float');  % Landmark 
hdr.rdb.last_slice_num = fread(id, 1, 'int');   
hdr.rdb.asset_slice_R = fread(id, 1, 'float');  % Slice reduction factor 
hdr.rdb.asset_slabwrap = fread(id, 1, 'float');  
hdr.rdb.dwnav_coeff = fread(id, 1, 'float');  % Coeff for amount of phase correction 
hdr.rdb.dwnav_cor = fread(id, 1, 'short');  % Navigator echo correction 
hdr.rdb.dwnav_view = fread(id, 1, 'short');  % Num of views of nav echoes 
hdr.rdb.dwnav_corecho = fread(id, 1, 'short');  % Num of nav echoes for actual correction 
hdr.rdb.dwnav_sview = fread(id, 1, 'short');  % Start view for phase correction process 
hdr.rdb.dwnav_eview = fread(id, 1, 'short');  % End view for phase correction process 
hdr.rdb.dwnav_sshot = fread(id, 1, 'short');  % Start shot for delta phase estimation in nav echoes 
hdr.rdb.dwnav_eshot = fread(id, 1, 'short');  % End shot for delta phase estimation in nav echoes 
hdr.rdb.win3d_type = fread(id, 1, 'short');  % 0 = Modified Hanning, 1 = modified Tukey 
hdr.rdb.win3d_apod = fread(id, 1, 'float');  % degree of apodization; 0.0 = boxcar, 1.0=hanning 
hdr.rdb.win3d_q = fread(id, 1, 'float');  % apodization at ends, 0.0 = max, 1.0 = boxcar 
hdr.rdb.ime_scic_enable = fread(id, 1, 'short');  % Surface Coil Intensity Correction: 1 if enabled 
hdr.rdb.clariview_type = fread(id, 1, 'short');  % Type of Clariview/Name of Filter 
hdr.rdb.ime_scic_edge = fread(id, 1, 'float');  % Edge paramaters for Enhanced Recon 
hdr.rdb.ime_scic_smooth = fread(id, 1, 'float');  % Smooth paramaters for Enhanced Recon 
hdr.rdb.ime_scic_focus = fread(id, 1, 'float');  % Focus paramaters for Enhanced Recon 
hdr.rdb.clariview_edge = fread(id, 1, 'float');  % Edge paramaters for clariview 
hdr.rdb.clariview_smooth = fread(id, 1, 'float');  % Smooth paramaters for clariview 
hdr.rdb.clariview_focus = fread(id, 1, 'float');  % Focus paramaters for clariview 
hdr.rdb.scic_reduction = fread(id, 1, 'float');  % Reduction paramater for SCIC 
hdr.rdb.scic_gauss = fread(id, 1, 'float');  % Gauss paramater for SCIC 
hdr.rdb.scic_threshold = fread(id, 1, 'float');  % Threshold paramater for SCIC 
hdr.rdb.ectricks_no_regions = fread(id, 1, 'long');  % Total no of regions acquired by PSD 
hdr.rdb.ectricks_input_regions = fread(id, 1, 'long');  % Total no of input regions for reordering 
hdr.rdb.psc_reuse = fread(id, 1, 'short');  % Header field for smart prescan 
hdr.rdb.left_blank = fread(id, 1, 'short');  
hdr.rdb.right_blank = fread(id, 1, 'short');  
hdr.rdb.acquire_type = fread(id, 1, 'short');  % Acquire type information from CV 
hdr.rdb.retro_control = fread(id, 1, 'short');  % Retrosective FSE phase correction control flag. (continued...)
hdr.rdb.etl = fread(id, 1, 'short');  % Added for Retrospective FSE phase correction. This (continued...)
hdr.rdb.pcref_start = fread(id, 1, 'short');  % 1st view to use for dynamic EPI phase correction. 
hdr.rdb.pcref_stop = fread(id, 1, 'short');  % Last view to use for dynamic EPI phase correction. 
hdr.rdb.ref_skip = fread(id, 1, 'short');  % Number of passes to skip for dynamic EPI phase correction. 
hdr.rdb.extra_frames_top = fread(id, 1, 'short');  % Number of extra frames at top of K-space 
hdr.rdb.extra_frames_bot = fread(id, 1, 'short');  % Number of extra frames at bottom of K-space 
hdr.rdb.multiphase_type = fread(id, 1, 'short');  % 0 = INTERLEAVED , 1 = SEQUENTIAL 
hdr.rdb.nphases = fread(id, 1, 'short');  % Number of phases in a multiphase scan 
hdr.rdb.pure = fread(id, 1, 'short');  % PURE flag from psd 
hdr.rdb.pure_scale = fread(id, 1, 'float');  % Recon scale factor ratio for cal scan 
hdr.rdb.off_data = fread(id, 1, 'int');  % Byte offset to start of raw data (i.e size of POOL_HEADER) 
hdr.rdb.off_per_pass = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_per_pass of POOL_HEADER 
hdr.rdb.off_unlock_raw = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_unlock_raw of POOL_HEADER 
hdr.rdb.off_data_acq_tab = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_data_acq_tab of POOL_HEADER 
hdr.rdb.off_nex_tab = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_nex_tab of POOL_HEADER 
hdr.rdb.off_nex_abort_tab = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_nex_abort_tab of POOL_HEADER 
hdr.rdb.off_tool = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_tool of POOL_HEADER 
hdr.rdb.off_exam = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_exam of POOL_HEADER 
hdr.rdb.off_series = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_series of POOL_HEADER 
hdr.rdb.off_image = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_image of POOL_HEADER 
hdr.rdb.off_ps = fread(id, 1, 'int');  % Byte offset to start of rdb_hdr_ps of POOL_HEADER 
hdr.rdb.off_spare_b = fread(id, 1, 'int');  % spare 
hdr.rdb.new_wnd_level_flag = fread(id, 1, 'int');  % New WW/WL algo enable/disable flag 
hdr.rdb.wnd_image_hist_area = fread(id, 1, 'int');  % Image Area % 
hdr.rdb.wnd_high_hist = fread(id, 1, 'float');  % Histogram Area Top 
hdr.rdb.wnd_lower_hist = fread(id, 1, 'float');  % Histogram Area Bottom 
hdr.rdb.pure_filter = fread(id, 1, 'short');  % PURE noise reduction on=1/off=0 
hdr.rdb.cfg_pure_filter = fread(id, 1, 'short');  % PURE cfg file value 
hdr.rdb.cfg_pure_fit_order = fread(id, 1, 'short');  % PURE cfg file value 
hdr.rdb.cfg_pure_kernelsize_z = fread(id, 1, 'short');  % PURE cfg file value 
hdr.rdb.cfg_pure_kernelsize_xy = fread(id, 1, 'short');  % PURE cfg file value 
hdr.rdb.cfg_pure_weight_radius = fread(id, 1, 'short');  % PURE cfg file value 
hdr.rdb.cfg_pure_intensity_scale = fread(id, 1, 'short');  % PURE cfg file value 
hdr.rdb.cfg_pure_noise_threshold = fread(id, 1, 'short');  % PURE cfg file value 
hdr.rdb.wienera = fread(id, 1, 'float');  % NB maintain alignment of floats 
hdr.rdb.wienerb = fread(id, 1, 'float');  
hdr.rdb.wienert2 = fread(id, 1, 'float');  
hdr.rdb.wieneresp = fread(id, 1, 'float');  
hdr.rdb.wiener = fread(id, 1, 'short');  
hdr.rdb.flipfilter = fread(id, 1, 'short');  
hdr.rdb.dbgrecon = fread(id, 1, 'short');  
hdr.rdb.ech2skip = fread(id, 1, 'short');  
hdr.rdb.tricks_type = fread(id, 1, 'int');  % 0 = Subtracted, 1 = Unsubtracted 
hdr.rdb.lcfiesta_phase = fread(id, 1, 'float');  % LC Fiesta 
hdr.rdb.lcfiesta = fread(id, 1, 'short');  % LC Fiesta 
hdr.rdb.herawflt = fread(id, 1, 'short');  % Half echo raw data filter 
hdr.rdb.herawflt_befnwin = fread(id, 1, 'short');  % Half echo raw data filter 
hdr.rdb.herawflt_befntran = fread(id, 1, 'short');  % Half echo raw data filter 
hdr.rdb.herawflt_befamp = fread(id, 1, 'float');  % Half echo raw data filter 
hdr.rdb.herawflt_hpfamp = fread(id, 1, 'float');  % Half echo raw data filter 
hdr.rdb.heover = fread(id, 1, 'short');  % Half echo over sampling 
hdr.rdb.pure_correction_threshold = fread(id, 1, 'short');  % PURE Correction threshold 
hdr.rdb.swiftenable = fread(id, 1, 'int');  % SWIFT enable/disable flag 
hdr.rdb.numslabs = fread(id, 1, 'short');  % Number of slabs to be used by TRICKS 
hdr.rdb.swiftcoilnos = fread(id, 1, 'ushort');  % Number of coils to SWIFT between 
hdr.rdb.ps_autoshim_status = fread(id, 1, 'int');   
hdr.rdb.dynaplan_numphases = fread(id, 1, 'int');  % Number of phases for Dynamic Plan 
hdr.rdb.medal_cfg = fread(id, 1, 'short');  % MEDAL configuration bitmask 
hdr.rdb.medal_nstack = fread(id, 1, 'short');  % MEDAL pixel stack size 
hdr.rdb.medal_echo_order = fread(id, 1, 'short');  % MEDAL in-phase, out-of-phase order
hdr.rdb.medal_kernel_up = fread(id, 1, 'short');  % MEDAL max adaptive region (continued...)
hdr.rdb.medal_kernel_down = fread(id, 1, 'short');  % MEDAL min adaptive region (continued...)
hdr.rdb.medal_kernel_smooth = fread(id, 1, 'short');  % MEDAL field smoothing (continued...)
hdr.rdb.medal_start = fread(id, 1, 'short');  % MEDAL 2D region growing (continued...)
hdr.rdb.medal_end = fread(id, 1, 'short');  % MEDAL 2D region growing (continued...)
hdr.rdb.rcideal = fread(id, 1, 'uint');  % Enable/Disable flag 
hdr.rdb.rcdixproc = fread(id, 1, 'uint');  % IDEAL image options and IDEAL control bits 
hdr.rdb.df = fread(id, 1, 'float');  % Delta Frequency between two species 
hdr.rdb.bw = fread(id, 1, 'float');  % Bandwidth in hz 
hdr.rdb.te1 = fread(id, 1, 'float');  % First Echo time in ms 
hdr.rdb.esp = fread(id, 1, 'float');  % Echo spacing in ms 
hdr.rdb.feextra = fread(id, 1, 'int');  % This will give the number of extra points 
hdr.rdb.vibrant = fread(id, 1, 'short');  % Set to 1 for VIBRANT scans 
hdr.rdb.asset_torso = fread(id, 1, 'short');  % Set to 1 for torso scans 
hdr.rdb.asset_alt_cal = fread(id, 1, 'int');  % Set to 1 to use reapodized cal 
hdr.rdb.kacq_uid = fread(id, 1, 'int');  % unique id for kacq_yz.txt files 
hdr.rdb.asset_localTx = fread(id, 1, 'int');  % Set to 1 for local Tx phased arrays 
hdr.rdb.threedscale = fread(id, 1, 'float');  % modified: aglatz; Use to scale 3D acqs to avoid overrange 
hdr.rdb.scanner_mode = fread(id, 1, 'short');  % 1=Product, 2=Research, 3=Service 
hdr.rdb.short_padding_1 = fread(id, 1, 'short');  
hdr.rdb.excess = fread(id, 182, 'short');  % free space for later expansion 
hdr.per_pass = fread(id, 4096, 'char');  % lumped type RDB_PER_PASS_TAB
hdr.unlock_raw = fread(id, 4096, 'char');  % lumped type RDB_PER_PASS_TAB
for islice=[1:1024]  % kluge for type RDB_SLICE_INFO_ENTRY
  hdr.data_acq_tab(islice).pass_number = fread(id, 1, 'short');  % which pass this slice is in
  hdr.data_acq_tab(islice).slice_in_pass = fread(id, 1, 'short');  % which slice in this pass
  hdr.data_acq_tab(islice).gw_point = fread(id, [3,3], 'float');  % corner points of image
end
hdr.nex_tab = fread(id, 2052, 'char');  % lumped type RDB_NEX_TYPE
hdr.nex_abort_tab = fread(id, 2052, 'char');  % lumped type RDB_NEX_TYPE
hdr.tool = fread(id, 2048, 'char');  % lumped type TOOLSDATA
