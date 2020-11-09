function [times dkangle sched tt1 tt2] = ffbm_findbmruntimes(experiment,mapsourcetype,year,rxNum)
%function [times dkangle sched tt1 tt2] = ffbm_findbmruntimes(experiment,mapsourcetype,year,rxNum)
%
% Function to find the beam map run start and stop times
% The function is handwritten, so maybe we should make
% a csv file instead, but this has worked so far.
% As more years are added, will need to add!
%
% Inputs:
%         experiment:    'bicep2','keck'
%         mapsourcetype: 'uberchopper','zechoppa' for bicep2
%         year:          '2012','2013','2014' for keck
%         rxNum:         0, 1, 2, 3, 4
%                        (only for keck, can be left empty for Bicep2)
%                        'all' for all Keck rxs
%
% CLW 20140515
% KSK 20140525 added 2014 BM runs

switch experiment
  
%%%%%%%% BICEP2 %%%%%%%%%%%%%%%%%%
  case 'bicep2'
    
    switch year
      case 2012
	mapsourcetype = 'uberchopper';
      case 2011
	mapsourcetype = 'zechoppa';
    end
        
    switch mapsourcetype
      case 'uberchopper'
	
	time={};
	
	time{1}='20121113_043509';%dk 90
	t1{1}='2012-nov-13:04:36:39';
	t2{1}='2012-nov-13:13:41:57';
	
	time{2}='20121113_134005';%dk 0
	t1{2}='2012-nov-13:13:42:00';
	t2{2}='2012-nov-13:22:47:24';
	
	time{3}='20121113_224436';%dk -90
	t1{3}='2012-nov-13:22:47:28';
	t2{3}='2012-nov-14:07:52:51';
	
	time{4}='20121114_075224';%dk -180
	t1{4}='2012-nov-14:07:52:56';
	t2{4}='2012-nov-14:16:58:47';
	
	time{5}='20121115_054553';%90
	t1{5}='2012-nov-15:05:47:25';
	t2{5}='2012-nov-15:14:52:30';
	
	time{6}='20121115_145049';%0
	t1{6}='2012-nov-15:14:52:33';
	t2{6}='2012-nov-15:23:57:56';
	
	time{7}='20121116_014434';%-90
	t1{7}='2012-nov-16:01:46:06';
	t2{7}='2012-nov-16:10:51:45';
	
	time{8}='20121116_104930';%-180
	t1{8}='2012-nov-16:10:51:49';
	t2{8}='2012-nov-16:19:57:41';
	
	time{9}='20121117_022104';%90
	t1{9}='2012-nov-17:02:22:36';
	t2{9}='2012-nov-17:11:28:13';
	
	time{10}='20121117_112600';%0
	t1{10}='2012-nov-17:11:28:17';
	t2{10}='2012-nov-17:20:33:40';
	
	time{11}='20121118_045302';%-90
	t1{11}='2012-nov-18:04:54:33';
	t2{11}='2012-nov-18:13:59:40';
	
	time{12}='20121118_135758';%-180
	t1{12}='2012-nov-18:13:59:45';
	t2{12}='2012-nov-18:23:05:36';
	
	dkangle=[90 0 -90 -180 90 0 -90 -180 90 0 -90 -180];
	sched=[26 26 26 26 27 27 27 27 28 28 28 28];
	
      case 'zechoppa'
	
	time={};
	
	time{1}='20111120_235000';%dk 90
	t1{1}='2011-nov-20:23:50';
	t2{1}='2011-nov-21:04:39';

	%20111121_045044
	jj=2;
	time{jj}='20111121_044726';%dk 0 
	t1{jj}='2011-nov-21:04:50';
	t2{jj}='2011-nov-21:09:28';
		
	%20111121_092757
	jj=3;
	time{jj}='20111121_092757';%dk-90  
	t1{jj}='2011-nov-21:09:30';
	t2{jj}='2011-nov-21:14:15';
	
	%20111121_141503
	jj=4;
	time{jj}='20111121_141503';%dk-180  
	t1{jj}='2011-nov-21:14:18';
	t2{jj}='2011-nov-21:19:05';

	%20111202_020734
	jj=5;
	time{jj}='20111202_020734';%dk  90
	t1{jj}='2011-dec-02:02:10';
	t2{jj}='2011-dec-02:06:56';

	%20111202_065506
	jj=6;
	time{jj}='20111202_065824';%dk 0
	t1{jj}='2011-dec-02:07:00';
	t2{jj}='2011-dec-02:11:45';

	%20111202_114530
	jj=7;
	time{jj}='20111202_114530';%dk -90
	t1{jj}='2011-dec-02:11:47';
	t2{jj}='2011-dec-02:16:34';

	%20111202_163236
	jj=8;
	time{jj}='20111202_163554';%dk -180
	t1{jj}='2011-dec-02:16:36';
	t2{jj}='2011-dec-02:21:23';

	%20111202_212839
	jj=9;
	time{jj}='20111202_212839';%dk 90
	t1{jj}='2011-dec-02:21:27';
	t2{jj}='2011-dec-03:02:14';
	
	%20111203_073234
	jj=10;
	time{jj}='20111203_073234';%dk 0
	t1{jj}='2011-dec-03:07:33';
	t2{jj}='2011-dec-03:12:22';
		
	%20111203_122006
	jj=11;
	time{jj}='20111203_122006';%dk -90
	t1{jj}='2011-dec-03:12:23';
	t2{jj}='2011-dec-03:17:10';

	%20111203_171030
	jj=12;
	time{jj}='20111203_171030';%dk  
	t1{jj}='2011-dec-03:17:12';
	t2{jj}='2011-dec-03:22:00';
	
	dkangle=[90 0 -90 -180 90 0 -90 -180 90 0 -90 -180];
	sched=[1 1 1 1 2 2 2 2 3 3 3 3];
    
    end
    
    times=time;
    tt1=t1;
    tt2=t2;


%%%%%%%% KECK %%%%%%%%%%%%%%%%%%
  case 'keck'
    
    switch year
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%
      case 2012
	%7_ffflat_dslmast_28a_000.sch, rx0,1,4
	time{1}='20120213_053604';
	t1{1}='2012-feb-13:05:31:28';
	t2{1}='2012-feb-13:13:12:12';
	dk(1)=-122;
	sc(1)=0;
	
	%7_ffflat_dslmast_28e_000.sch ,rx3,4,0
	t1{2}='2012-feb-14:21:51:03';
	t2{2}='2012-feb-15:05:31:44';
	time{2}='20120214_215537';
	dk(2)=-194;
	sc(2)=0;
	
	%7_ffflat_dslmast_28a_001.sch rx4,0,1
	jj=3;
	t1{jj}='2012-feb-15:05:31:44'; 
	t2{jj}='2012-feb-15:13:12:22';
	time{jj}='20120215_053015';
	dk(jj)=-122;
	sc(jj)=1;
	
	%7_ffflat_dslmast_28b_001.sch rx0,1,2
	jj=4;
	t1{jj}='2012-feb-15:13:12:23'; 
	t2{jj}='2012-feb-15:20:52:51'; 
	time{jj}='20120215_131254';
	dk(jj)=-50;
	sc(jj)=1;
	
	%7_ffflat_dslmast_28c_001.sch  rx1,2,3
	jj=5;
	t1{jj}='2012-feb-15:20:52:52'; 
	t2{jj}='2012-feb-16:04:33:47'; 
	time{jj}='20120215_205726';
	dk(jj)=22;
	sc(jj)=1;
	
	%7_ffflat_dslmast_28d_001.sch  rx2,3,4
	jj=6;
	t1{jj}='2012-feb-16:14:51:22'; 
	t2{jj}='2012-feb-16:22:30:23'; 
	time{jj}='20120216_145114';
	dk(jj)=-266;%+94
	sc(jj)=1;
	
	%7_ffflat_dslmast_28e_001.sch rx3,4,0
	jj=7;
	t1{jj}='2012-feb-16:22:30:24'; 
	t2{jj}='2012-feb-17:06:11:09'; 
	time{jj}='20120216_223055';
	dk(jj)=-194;
	sc(jj)=1;
	
	%7_ffflat_dslmast_28a_002.sch rx4,0,1
	jj=8;
	t1{jj}='2012-feb-17:06:11:10'; 
	t2{jj}='2012-feb-17:13:51:43'; 
	time{jj}='20120217_061543';
	dk(jj)=-122;
	sc(jj)=2;
	
	%7_ffflat_dslmast_28b_002.sch rx0,1,2
	jj=9;
	t1{jj}='2012-feb-17:13:51:43'; 
	t2{jj}='2012-feb-17:21:32:16'; 
	time{jj}='20120217_135214';
	dk(jj)=-50;
	sc(jj)=2;
	
	%7_ffflat_dslmast_28c_002.sch rx1,2,3
	jj=10;
	t1{jj}='2012-feb-18:03:26:46'; 
	t2{jj}='2012-feb-18:11:07:36'; 
	time{jj}='20120218_033120';
	dk(jj)=22;
	sc(jj)=2;
	
	%7_ffflat_dslmast_28d_002.sch rx2,3,4
	jj=11;
	t1{jj}='2012-feb-18:11:07:37'; 
	t2{jj}='2012-feb-18:18:50:14'; 
	time{jj}='20120218_111211';
	dk(jj)=-266;
	sc(jj)=2;
	
	%7_ffflat_dslmast_28e_002.sch rx3,4,0
	jj=12;
	t1{jj}='2012-feb-18:18:50:15'; 
	t2{jj}='2012-feb-19:02:30:42'; 
	time{jj}='20120218_185448';
	dk(jj)=-194;
	sc(jj)=2;
	
	%7_ffflat_dslmast_28c_000.sch rx1,2,3
	jj=13;
	t1{jj}='2012-feb-19:02:30:43'; 
	t2{jj}='2012-feb-19:10:11:52'; 
	time{jj}='20120219_023517';
	dk(jj)=22;
	sc(jj)=0;
	
	%7_ffflat_dslmast_28d_003.sch rx2,3,4
	jj=14;
	t1{jj}='2012-feb-21:03:49:41'; 
	t2{jj}='2012-feb-21:11:29:54'; 
	time{jj}='20120221_035415';
	dk(jj)=-266;
	sc(jj)=0;
	
	%7_ffflat_dslmast_28i_000.sch rx3,4
	jj=15;
	t1{jj}='2012-feb-21:11:29:55'; 
	t2{jj}='2012-feb-21:19:09:59'; 
	time{jj}='20120221_113429';
	dk(jj)=-230;
	sc(jj)=1;
	
	%7_ffflat_dslmast_28j_000.sch rx4,0
	jj=16;
	t1{jj}='2012-feb-21:19:10:00'; 
	t2{jj}='2012-feb-22:02:50:19'; 
	time{jj}='20120221_191434';
	dk(jj)=-158;
	sc(jj)=1;
	
	%7_ffflat_dslmast_28h_000.sch rx2,3
	jj=17;
	t1{jj}='2012-feb-22:08:25:04'; 
	t2{jj}='2012-feb-22:16:05:42';
	time{jj}='20120222_082938';
	dk(jj)=-302;
	sc(jj)=1;
	
	%7_ffflat_dslmast_28f_000.sch rx0,1
	jj=18;
	t1{jj}='2012-feb-22:16:05:42'; 
	t2{jj}='2012-feb-22:23:47:46'; 
	time{jj}='20120222_161015';
	dk(jj)=-86;
	sc(jj)=1;
	
	%7_ffflat_dslmast_28g_000.sch rx1 2
	jj=19;
	t1{jj}='2012-feb-22:23:47:47'; 
	t2{jj}='2012-feb-23:07:28:25'; 
	time{jj}='20120222_235221';
	dk(jj)=-14;
	sc(jj)=1;
	
	%7_ffflat_dslmast_28f_004.sch rx0,1
	jj=20;
	t1{jj}='2012-feb-24:11:26:20'; 
	t2{jj}='2012-feb-24:19:06:30'; 
	time{jj}='20120224_113054';
	dk(jj)=-86;
	sc(jj)=2;
	
	%7_ffflat_dslmast_28g_004.sch rx1,2
	jj=21;
	t1{jj}='2012-feb-24:19:06:31'; 
	t2{jj}='2012-feb-25:02:47:23'; 
	time{jj}='20120224_191105';
	dk(jj)=-14;
	sc(jj)=2;
	
	%7_ffflat_dslmast_28h_004.sch rx2,3
	jj=22;
	t1{jj}='2012-feb-25:08:18:20'; 
	t2{jj}='2012-feb-25:15:59:44'; 
	time{jj}='20120225_082234';
	dk(jj)=-302;
	sc(jj)=2;
	
	%7_ffflat_dslmast_28j_005.sch rx4,0
	jj=23;
	t1{jj}='2012-feb-26:04:01:00'; 
	t2{jj}='2012-feb-26:11:39:18'; 
	time{jj}='20120226_040412';
	dk(jj)=-158;
	sc(jj)=2;
	
	%7_ffflat_dslmast_28i_005.sch rx3,4
	jj=24;
	t1{jj}='2012-feb-27:00:46:52'; 
	t2{jj}='2012-feb-27:08:27:37'; 
	time{jj}='20120227_005126';
	dk(jj)=-230;
	sc(jj)=2;
	
	%7_ffflat_dslmast_28h_005.sch rx2,3
	jj=25;
	t1{jj}='2012-feb-27:08:27:38'; 
	t2{jj}='2012-feb-27:16:08:17';
	time{jj}='20120227_083213';
	dk(jj)=-302;
	sc(jj)=0;
	
	switch rxNum
	  case 0
	    index=[1 2 3 4 7 8 9 12 16 18 20 23];%index([8])=[];
	  case 1
	    index=[1 3 4 5 8 9 10 13 18 19 20 21];
	    %index=[1:length(time)];
	  case 2
	    index=[4 5 6 9 10 11 13 14 17 19 21 22 25];index(12)=[];
	    %map 22 has a lot of stripes
	  case 3
	    index=[2 5 6 7 10 11 12 13 14 15 17 22 24 25];index(12)=[];
	  case 4
	    index=[1 2 3 6 7 8 11 12 14 15 16 23 24];%index(8)=[];
	  case 'all'
	    index=[1:length(time)];%index(12)=[];
	end
	
	times={};
	tt1={};
	tt2={};
	dkangle=[];
	sched=[];
	
	for ii=1:length(index)
	  times=[times time{index(ii)}];
	  tt1=[tt1 t1{index(ii)}];
	  tt2=[tt2 t2{index(ii)}];
	  dkangle=[dkangle dk(index(ii))];
	  sched=[sched sc(index(ii))];
	end
	
      %%%%%%%%%%%%%%%%%%%%%%%%%%	
      case 2013
	
	time{1}='20130212_005232';%30b_000, rx1, rx2
	t1{1}='2013-feb-12:00:48:12';
	t2{1}='2013-feb-12:08:28:15';
	dk(1)=-50;
	sc(1)=0;
	
	time{2}='20130212_083241';%30c_000, rx1,2,3
	t1{2}='2013-feb-12:08:28:20';
	t2{2}='2013-feb-12:16:08:30';
	dk(2)=22;
	sc(2)=0;
	
	time{3}='20130213_014347';%30d_000, rx2,3,4
	t1{3}='2013-feb-13:01:39:26';
	t2{3}='2013-feb-13:09:19:40';
	dk(3)=-266;
	sc(3)=0;
	
	time{4}='20130213_092408';%30e_000, rx3,4
	t1{4}='2013-feb-13:09:19:47';
	t2{4}='2013-feb-13:16:59:34';
	dk(4)=-194;
	sc(4)=0;
	
	time{5}='20130219_122110';%31j_000, rx4
	t1{5}='2013-feb-19:12:16:50';
	t2{5}='2013-feb-19:19:57:00';
	dk(5)=-158;
	sc(5)=0;
	
	time{6}='20130219_212859';%31f_000, rx1
	t1{6}='2013-feb-19:21:24:44';
	t2{6}='2013-feb-20:05:05:12';
	dk(6)=-86;
	sc(6)=0;
	
	time{7}='20130220_211148';%31i_000, rx3,4
	t1{7}='2013-feb-20:21:07:27';
	t2{7}='2013-feb-21:04:46:47';
	dk(7)=-230;
	sc(7)=0;
	
	time{8}='20130221_045727';%31h_000,rx2,3
	t1{8}='2013-feb-21:04:53:06';
	t2{8}='2013-feb-21:12:32:52';
	dk(8)=-302;
	sc(8)=0;
	
	time{9}='20130221_123717';%31e_000,rx3,4
	t1{9}='2013-feb-21:12:32:56';
	t2{9}='2013-feb-21:20:13:42';
	dk(9)=-194;
	sc(9)=1;
	
	time{10}='20130221_201807';%31g_000,rx1,2
	t1{10}='2013-feb-21:20:14:16';
	t2{10}='2013-feb-22:03:55:16';
	dk(10)=-14;
	sc(10)=0;
	
	time{11}='20130224_031531';%31d_001,rx2,3,4
	t1{11}='2013-feb-24:03:11:16';
	t2{11}='2013-feb-24:10:50:29';
	dk(11)=-266;%+94
	sc(11)=1;
	
	time{12}='20130224_110215';%31e_001,rx0,3,4
	t1{12}='2013-feb-24:10:50:32';
	t2{12}='2013-feb-24:18:30:43';
	dk(12)=-194;
	sc(12)=2;
	
	time{13}='20130224_183116';%31h_001,rx2,3
	t1{13}='2013-feb-24:18:30:47';
	t2{13}='2013-feb-25:02:11:30';
	dk(13)=-302;
	sc(13)=1;
	
	time{14}='20130226_050440';%31g_001,rx1,2
	t1{14}='2013-feb-26:05:00:18';
	t2{14}='2013-feb-26:12:39:58';
	dk(14)=-14;
	sc(14)=1;
	
	time{15}='20130226_124424';%31c_001,rx1,2,3
	t1{15}='2013-feb-26:12:40:33';
	t2{15}='2013-feb-26:20:19:32';
	dk(15)=22;
	sc(15)=1;
	
	time{16}='20130301_142525';%31a_003,rx1,rx4, incomplete files
	t1{16}='2013-mar-01:14:26:00';
	t2{16}='2013-mar-01:22:04:36';
	dk(16)=-122;
	sc(16)=0;
	
	time{17}='20130301_220903';%31d_003,rx2,3,4,incomplete files
	t1{17}='2013-mar-01:22:05:11';
	t2{17}='2013-mar-02:05:45:47';
	dk(17)=-266;%+94
	sc(17)=2;
	
	time{18}='20130302_055012';%31i_003,rx3,4
	t1{18}='2013-mar-02:05:45:50';
	t2{18}='2013-mar-02:13:26:07';
	dk(18)=-230;
	sc(18)=1;
	
	time{19}='20130303_120516';%31h_003,rx2,3
	t1{19}='2013-mar-03:12:01:25';
	t2{19}='2013-mar-03:19:41:47';
	dk(19)=-302;
	sc(19)=2;
	
	time{20}='20130303_194955';%31b_003,rx1,2 
	t1{20}='2013-mar-03:19:42:22';
	t2{20}='2013-mar-04:03:24:51';
	dk(20)=-50;
	sc(20)=1;
	
	time{21}='20130304_032916';%31c_003, rx1,2,3,missing files
	t1{21}='2013-mar-04:03:24:54';
	t2{21}='2013-mar-04:11:05:43';
	dk(21)=22;
	sc(21)=2;
	
	time{22}='20130305_200933';%31j_004, rx4,need to run this one
	t1{22}='2013-mar-05:20:05:11';
	t2{22}='2013-mar-06:03:45:19';
	dk(22)=-158;
	sc(22)=1;
	
	time{23}='20130306_095508';%31a_004,rx1,4
	t1{23}='2013-mar-06:09:50:51';
	t2{23}='2013-mar-06:17:31:19';
	dk(23)=-122;
	sc(23)=1;
	
	time{24}='20130307_064412';%31a_005,rx0,1,4
	t1{24}='2013-mar-07:06:39:50';
	t2{24}='2013-mar-07:14:19:22';
	dk(24)=-122;
	sc(24)=2;
	
	time{25}='20130307_142347';%31f_005,rx0,1
	t1{25}='2013-mar-07:14:19:25';
	t2{25}='2013-mar-07:21:59:02';
	dk(25)=-86;
	sc(25)=1;
	
	time{26}='20130308_062140';%31j_005,rx0,4
	t1{26}='2013-mar-08:06:17:19';
	t2{26}='2013-mar-08:13:57:55';
	dk(26)=-158;
	sc(26)=2;
	
	switch rxNum
	  case 0
	    index=[12 24 25 26];
	  case 1
	    index=[1 2 6 10 14 15 16 20 21 23 24 25];
	  case 2
	    index=[1 2 3 8 10 11 13 14 15 17 19 20 21];
	  case 3
	    index=[2 3 4 7 8 9 11 12 13 15 17 18 19 21];
	  case 4
	    index=[3 4 5 7 9 11 12 16 17 18 22 23 24 26];
	  case 'all'
	    index=[1:length(time)];%
	end
	
	times={};
	tt1={};
	tt2={};
	dkangle=[];
	sched=[];
	
	for ii=1:length(index)
	  times=[times time{index(ii)}];
	  tt1=[tt1 t1{index(ii)}];
	  tt2=[tt2 t2{index(ii)}];
	  dkangle=[dkangle dk(index(ii))];
	  sched=[sched sc(index(ii))];
	end

      %%%%%%%%%%%%%%%%%%%%%%%%%%
      case 2014
      
	%7_ffflat_dslmast_31a_003.sch
	jj=1;
	time{jj}='20140212_060550';
	t1{jj}='2014-feb-12:06:05:10';
	t2{jj}='2014-feb-12:13:44:40';
	dk(jj)=-122;
	sc(jj)=1;
	
	%7_ffflat_dslmast_31b_005.sch
	jj=2;
	time{jj}='20140212_204427';
	t1{jj}='2014-feb-12:20:43:50';
	t2{jj}='2014-feb-13:04:23:40';
	dk(jj)=-50;
	sc(jj)=1;
				

	%7_ffflat_dslmast_31c_006.sch
	jj=3;
	time{jj}='20140213_100548';
	t1{jj}='2014-feb-13:10:05:15';
	t2{jj}='2014-feb-13:17:45:00';
	dk(jj)=22;
	sc(jj)=1;
	
	%7_ffflat_dslmast_31d_006.sch
	jj=4;
	time{jj}='20140213_174943';
	t1{jj}='2014-feb-13:17:45:13';
	t2{jj}='2014-feb-14:01:26:50';
	dk(jj)=-266;
	sc(jj)=1;
	
	%choppe freq is 14Hz
	%7_ffflat_dslmast_31e_006.sch
	jj=5;
	time{jj}='20140214_022039';
	t1{jj}='2014-feb-14:02:16:10';
	t2{jj}='2014-feb-14:09:56:00';
	dk(jj)=-194;
	sc(jj)=1;
	
	%7_ffflat_dslmast_31f_006.sch
	jj=6;
	time{jj}='20140214_100044';
	t1{jj}='2014-feb-14:09:56:15';
	t2{jj}='2014-feb-14:17:36:50';
	dk(jj)=-86;
	sc(jj)=1;

	%7_ffflat_dslmast_31g_006.sch
	jj=7;
	time{jj}='20140215_024810';
	t1{jj}='2014-feb-15:02:43:45';
	t2{jj}='2014-feb-15:10:24:41';
	dk(jj)=-14;
	sc(jj)=1;
	
	%7_ffflat_dslmast_31h_006.sch
	jj=8;
	time{jj}='20140215_102917';
	t1{jj}='2014-feb-15:10:24:50';
	t2{jj}='2014-feb-15:18:06:20';
	dk(jj)=-302;
	sc(jj)=1;
	
	%7_ffflat_dslmast_31i_006.sch
	jj=9;
	time{jj}='20140215_181103';
	t1{jj}='2014-feb-15:18:06:30';
	t2{jj}='2014-feb-16:01:46:50';
	dk(jj)=-230;
	sc(jj)=1;
	
	%7_ffflat_dslmast_31j_006.sch
	jj=10;
	time{jj}='20140216_015133';
	t1{jj}='2014-feb-16:01:47:00';
	t2{jj}='2014-feb-16:09:27:26';
	dk(jj)=-158;
	sc(jj)=1;
	
	%7_ffflat_dslmast_31b_007.sch
	jj=11;
	time{jj}='20140216_153718';
	t1{jj}='2014-feb-16:15:32:43';
	t2{jj}='2014-feb-16:23:14:11';
	dk(jj)=-50;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31a_007.sch 
	% might need to be repeated
	%aborted
	%jj=12;
	%time{jj}='tmp';
	%t1{jj}='2014-feb-16:23:14:12';
	%t2{jj}='2014-feb-17:00:13:56';
	%dk(jj)=500;
	%sc(jj)=2;
	
	%7_ffflat_dslmast_31a_007.sch %repeated
	jj=12;
	time{jj}='20140217_001454';
	t1{jj}='2014-feb-17:00:15:00';
	t2{jj}='2014-feb-17:07:54:32';
	dk(jj)=-122;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31C_007.sch %repeated
	jj=13;
	time{jj}='20140217_075908';
	t1{jj}='2014-feb-17:07:54:40';
	t2{jj}='2014-feb-17:15:35:40';
	dk(jj)=22;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31d_007.sch %repeated
	jj=14;
	time{jj}='20140217_154023';
	t1{jj}='2014-feb-17:15:35:50';
	t2{jj}='2014-feb-17:23:17:30';
	dk(jj)=-266;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31e_007.sch %repeated
	jj=15;
	time{jj}='20140218_080933';
	t1{jj}='2014-feb-18:08:08:55';
	t2{jj}='2014-feb-18:15:48:28';
	dk(jj)=-194;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31f_007.sch %repeated
	jj=16;
	time{jj}='20140218_154747';
	t1{jj}='2014-feb-18:15:48:28';
	t2{jj}='2014-feb-18:23:29:21';
	dk(jj)=-86;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31g_007.sch %repeated
	jj=17;
	time{jj}='20140218_233005';
	t1{jj}='2014-feb-18:23:29:30';
	t2{jj}='2014-feb-19:07:09:50';
	dk(jj)=-14;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31h_007.sch %repeated
	jj=18;
	time{jj}='20140219_071426';
	t1{jj}='2014-feb-19:07:09:51';
	t2{jj}='2014-feb-19:14:51:36';
	dk(jj)=-302;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31i_007.sch
	jj=19;
	time{jj}='20140219_203955';
	t1{jj}='2014-feb-19:20:39:20';
	t2{jj}='2014-feb-20:04:19:10';
	dk(jj)=-230;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31j_007.sch
	jj=20;
	time{jj}='20140220_041958';
	t1{jj}='2014-feb-20:04:19:20';
	t2{jj}='2014-feb-20:11:59:30';
	dk(jj)=-158;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31b_008.sch
	jj=21;
	time{jj}='20140220_120023';
	t1{jj}='2014-feb-20:11:59:45';
	t2{jj}='2014-feb-20:19:40:30';
	dk(jj)=-50;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31a_008.sch
	jj=22;
	time{jj}='20140220_194118';
	t1{jj}='2014-feb-20:19:40:40';
	t2{jj}='2014-feb-21:03:21:20';
	dk(jj)=-122;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31c_008.sch
	jj=23;
	time{jj}='20140221_092016';
	t1{jj}='2014-feb-21:09:19:40';
	t2{jj}='2014-feb-21:17:00:00';
	dk(jj)=22;
	sc(jj)=2;

	%7_ffflat_dslmast_31d_008.sch
	jj=24;
	time{jj}='20140221_170444';
	t1{jj}='2014-feb-21:17:00:15';
	t2{jj}='2014-feb-22:00:41:50';
	dk(jj)=-266;
	sc(jj)=2;

	%7_ffflat_dslmast_31b_009.sch
	jj=25;
	time{jj}='20140222_024932';
	t1{jj}='2014-feb-22:02:45:00';
	t2{jj}='2014-feb-22:10:26:05';
	dk(jj)=-50;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31a_009.sch
	jj=26;
	time{jj}='20140222_102404';
	t1{jj}='2014-feb-22:10:26:05';
	t2{jj}='2014-feb-22:18:06:40';
	dk(jj)=-122;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31c_009.sch
	jj=27;
	time{jj}='20140222_235122';
	t1{jj}='2014-feb-22:23:50:50';
	t2{jj}='2014-feb-23:07:31:35';
	dk(jj)=22;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31d_009.sch
	jj=28;
	time{jj}='20140223_073616';
	t1{jj}='2014-feb-23:07:31:45';
	t2{jj}='2014-feb-23:15:13:20';
	dk(jj)=-266;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31e_009.sch
	jj=29;
	time{jj}='20140223_151409';
	t1{jj}='2014-feb-23:15:13:30';
	t2{jj}='2014-feb-23:22:54:00';
	dk(jj)=-194;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31f_009.sch
	jj=30;
	time{jj}='20140223_225449';
	t1{jj}='2014-feb-23:22:54:10';
	t2{jj}='2014-feb-24:06:35:10';
	dk(jj)=-86;
	sc(jj)=2;

	%7_ffflat_dslmast_31g_009.sch
	jj=31;
	time{jj}='20140224_115627';
	t1{jj}='2014-feb-24:11:55:55';
	t2{jj}='2014-feb-24:19:35:40';
	dk(jj)=-14;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31h_009.sch
	jj=32;
	time{jj}='20140224_194022';
	t1{jj}='2014-feb-24:19:35:47';
	t2{jj}='2014-feb-25:03:17:25';
	dk(jj)=-302;
	sc(jj)=2;

	%7_ffflat_dslmast_31i_009.sch
	jj=33;
	time{jj}='20140225_031815';
	t1{jj}='2014-feb-25:03:17:40';
	t2{jj}='2014-feb-25:10:57:50';
	dk(jj)=-230;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31j_009.sch
	jj=34;
	time{jj}='20140225_110234';
	t1{jj}='2014-feb-25:10:58:00';
	t2{jj}='2014-feb-25:18:38:40';
	dk(jj)=-158;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31a_010.sch
	jj=35;
	time{jj}='20140226_031317';
	t1{jj}='2014-feb-26:03:09:00';
	t2{jj}='2014-feb-26:10:48:40';
	dk(jj)=-122;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31b_010.sch
	jj=36;
	time{jj}='20140226_105325';
	t1{jj}='2014-feb-26:10:49:00';
	t2{jj}='2014-feb-26:18:29:00';
	dk(jj)=-50;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31c_010.sch
	jj=37;
	time{jj}='20140226_183344';
	t1{jj}='2014-feb-26:18:29:15';
	t2{jj}='2014-feb-27:02:09:50';
	dk(jj)=22;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31d_010.sch
	jj=38;
	time{jj}='20140227_021434';
	t1{jj}='2014-feb-27:02:10:05';
	t2{jj}='2014-feb-27:09:51:35';
	dk(jj)=-266;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31e_010.sch
	jj=39;
	time{jj}='20140227_193606';
	t1{jj}='2014-feb-27:19:31:45';
	t2{jj}='2014-feb-28:03:11:25';
	dk(jj)=-194;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31f_010.sch
	jj=40;
	time{jj}='20140228_031217';
	t1{jj}='2014-feb-28:03:11:40';
	t2{jj}='2014-feb-28:10:52:25';
	dk(jj)=-86;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31g_010.sch
	jj=41;
	time{jj}='20140228_105707';
	t1{jj}='2014-feb-28:10:52:40';
	t2{jj}='2014-feb-28:18:32:30';
	dk(jj)=-14;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31h_010.sch
	jj=42;
	time{jj}='20140228_183710';
	t1{jj}='2014-feb-28:18:32:45';
	t2{jj}='2014-mar-01:02:14:10';
	dk(jj)=-302;
	sc(jj)=2;

	%7_ffflat_dslmast_31i_010.sch
	jj=43;
	time{jj}='20140301_073616';
	t1{jj}='2014-mar-01:07:35:45';
	t2{jj}='2014-mar-01:15:15:25';
	dk(jj)=-230;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31j_010.sch
	jj=44;
	time{jj}='20140301_151619';
	t1{jj}='2014-mar-01:15:15:45';
	t2{jj}='2014-mar-01:22:55:30';
	dk(jj)=-158;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31h_011.sch
	jj=45;
	time{jj}='20140301_225622';
	t1{jj}='2014-mar-01:22:55:45';
	t2{jj}='2014-mar-02:06:36:35';
	dk(jj)=-302;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31g_011.sch
	jj=46;
	time{jj}='20140302_064502';
	t1{jj}='2014-mar-02:06:36:50';
	t2{jj}='2014-mar-02:14:19:20';
	dk(jj)=-14;
	sc(jj)=2;
	
	%7_ffflat_dslmast_31j_011.sch
	jj=47;
	time{jj}='20140302_231445';
	t1{jj}='2014-mar-02:23:10:15';
	t2{jj}='2014-mar-03:06:51:20';
	dk(jj)=-158;
	sc(jj)=2;

	%7_ffflat_dslmast_31a_013.sch
	jj=48;
	time{jj}='20140304_020019';
	t1{jj}='2014-mar-04:01:59:41';
	t2{jj}='2014-mar-04:09:39:21';
	dk(jj)=-122;
	sc(jj)=2;

	%7_ffflat_dslmast_31b_013.sch
	jj=49;
	time{jj}='20140304_094006';
	t1{jj}='2014-mar-04:09:39:22';
	t2{jj}='2014-mar-04:17:19:48';
	dk(jj)=-50;
	sc(jj)=2;

	%7_ffflat_dslmast_31c_013.sch
	jj=50;
	time{jj}='20140304_172032';
	t1{jj}='2014-mar-04:17:19:49';
	t2{jj}='2014-mar-05:01:00:39';
	dk(jj)=22;
	sc(jj)=2;

	%7_ffflat_dslmast_31d_013.sch
	jj=51;
	time{jj}='20140305_010515';
	t1{jj}='2014-mar-05:01:00:40';
	t2{jj}='2014-mar-05:08:42:29';
	dk(jj)=-266;
	sc(jj)=2;

	%7_ffflat_dslmast_31h_013.sch
	jj=52;
	time{jj}='20140305_141953';
	t1{jj}='2014-mar-05:14:15:22';
	t2{jj}='2014-mar-05:21:57:12';
	dk(jj)=-302;
	sc(jj)=2;

	%7_ffflat_dslmast_31g_013.sch
	jj=53;
	time{jj}='20140305_220148';
	t1{jj}='2014-mar-05:21:57:13';
	t2{jj}='2014-mar-06:05:38:59';
	dk(jj)=-14;
	sc(jj)=2;

	%7_ffflat_dslmast_31e_013.sch
	jj=54;
	time{jj}='20140306_053943';
	t1{jj}='2014-mar-06:05:39:00';
	t2{jj}='2014-mar-06:13:19:54';
	dk(jj)=-194;
	sc(jj)=2;

	%7_ffflat_dslmast_31f_013.sch
	jj=55;
	time{jj}='20140306_132038';
	t1{jj}='2014-mar-06:13:19:55';
	t2{jj}='2014-mar-06:21:00:14';
	dk(jj)=-86;
	sc(jj)=2;

	%7_ffflat_dslmast_31i_013.sch
	jj=56;
	time{jj}='20140307_021320';
	t1{jj}='2014-mar-07:02:12:42';
	t2{jj}='2014-mar-07:09:53:14';
	dk(jj)=-230;
	sc(jj)=2;

	%7_ffflat_dslmast_31j_013.sch
	jj=57;
	time{jj}='20140307_095358';
	t1{jj}='2014-mar-07:09:53:15';
	t2{jj}='2014-mar-07:17:33:40';
	dk(jj)=-158;
	sc(jj)=2;
	
	% Changed from back to front at beam map 24
	% Back: choose positions 1:5
	% Front: choose positions 1, 5:10
	switch rxNum
	  case 0
	    index = [01 02 05 06 10 11 12 15 16 20 21 22 ...
		     25 27 28 29 31 32 33 ...
		     36 37 38 39 41 42 43 ...
		     45 46 49 50 51 52 53 54 56];
	  case 1
	    index = [01 02 03 06 07 11 12 13 15 16 17 21 22 23 ...
		     26 27 28 29 32 33 34 ...
		     35 37 38 39 42 43 44 ...
		     45 47 48 50 51 52 54 56 57];
	  case 2
	    index = [02 03 04 07 08 11 13 14 17 18 21 23 24 ...
  		     25 26 28 29 30 33 34 ...
		     35 36 38 39 40 43 44 ...
		     47 48 49 51 54 55 56 57];
	  case 3
	    index = [03 04 05 08 09 13 14 18 19 23 24 ...
		     25 26 27 29 30 31 34 ...
		     35 36 37 39 40 41 44 ...
		     46 47 48 49 50 53 54 55 57];
	  case 4
	    index = [01 04 05 09 10 12 14 15 19 20 22 24 ...
		     25 26 27 28 30 31 32 ...
		     35 36 37 38 40 41 42 ...
		     45 46 48 49 50 51 52 53 55];
	  case 'all'
	    index = [1:length(time)];
	end
	
	times={};
	tt1={};
	tt2={};
	dkangle=[];
	sched=[];
	
	for ii=1:length(index)
	  times=[times time{index(ii)}];
	  tt1=[tt1 t1{index(ii)}];
	  tt2=[tt2 t2{index(ii)}];
	  dkangle=[dkangle dk(index(ii))];
	  sched=[sched sc(index(ii))];
	end
	
    end  % Year
    
end  % BICEP2 or KECK


