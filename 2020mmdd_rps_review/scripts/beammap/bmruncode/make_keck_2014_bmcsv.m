function make_keck_2014_bmccsv()
% Script is now obsolete. 
% Please use make_bmcsv which is more generic and yields the same result with 10x less work
% Script to generate a nice archival .csv for the 2014 Keck beam maps

disp('Script is now obsolete.')
disp(' Please use make_bmcsv with appropriate exp and year. It is more general  and yields the same result with 10x less work')
return  

% Don't overwrite a different experiment!
if ~strcmp(get_experiment_name(),'keck');
  error('Get into the right folder for Keck!')
end

savefile = 'beammaps/bmrunlist_2014.csv';

bm = bm_info();


for i = 1:length(bm)
  p.number(i) = bm(i).number;
  p.t1{i} = bm(i).t1;
  p.t2{i} = bm(i).t2;
  p.dk(i) = bm(i).dk;
  p.sched{i} = bm(i).sched;
  p.filename{i} = bm(i).filename;
  p.notes{i} = bm(i).notes;
  
  if i < 25
    p.mirror{i} = 'back';
  else
    p.mirror{i} = 'front';
  end
  
end

k.comments{1} = ['# BEAM MAPPING RUNS: KECK 2014'];
k.fields = {'number','t1','t2','dk','mirror','sched','filename','notes'};
k.formats = ...
    {'integer','string','string','double','string','string','string','string'};
k.units = ...
    {'unitless','time','time','deg','unitless','unitless','time','unitless'};


ParameterWrite(savefile,p,k);


return

function bm = bm_info()

bm(1).number = 1;
bm(1).t1 = '2014-feb-12:06:05:07';
bm(1).t2 = '2014-feb-12:13:44:45';
bm(1).dk = -122;
bm(1).sched = '7_ffflat_dslmast_31a_003.sch';
bm(1).filename = '20140212_060942';
bm(1).notes = '';

bm(2).number = 2;
bm(2).t1 = '2014-feb-12:20:43:44';
bm(2).t2 = '2014-feb-13:04:23:44';
bm(2).dk = -50;
bm(2).sched = '7_ffflat_dslmast_31b_005.sch';
bm(2).filename = '20140212_204819';
bm(2).notes = '';

bm(3).number = 3;
bm(3).t1 = '2014-feb-13:10:05:10';
bm(3).t2 = '2014-feb-13:17:45:07';
bm(3).dk = 22;
bm(3).sched = '7_ffflat_dslmast_31c_006.sch';
bm(3).filename = '20140213_100941';
bm(3).notes = '';

bm(4).number = 4;
bm(4).t1 = '2014-feb-13:17:45:08';
bm(4).t2 = '2014-feb-14:01:26:52';
bm(4).dk = 94;
bm(4).sched = '7_ffflat_dslmast_31d_006.sch';
bm(4).filename = '20140213_174943';
bm(4).notes = '';

bm(5).number = 5;
bm(5).t1 = '2014-feb-14:02:16:04';
bm(5).t2 = '2014-feb-14:09:56:08';
bm(5).dk = 166;
bm(5).sched = '7_ffflat_dslmast_31e_006.sch';
bm(5).filename = '20140214_022039';
bm(5).notes = 'Check Rx0';

bm(6).number = 6;
bm(6).t1 = '2014-feb-14:09:56:09';
bm(6).t2 = '2014-feb-14:17:36:55';
bm(6).dk = -86;
bm(6).sched = '7_ffflat_dslmast_31f_006.sch';
bm(6).filename = '20140214_100044';
bm(6).notes = 'Check Rx0';

bm(7).number = 7;
bm(7).t1 = '2014-feb-15:02:43:41';
bm(7).t2 = '2014-feb-15:10:24:41';
bm(7).dk = -14;
bm(7).sched = '7_ffflat_dslmast_31g_006.sch';
bm(7).filename = '20140215_024810';
bm(7).notes = '';

bm(8).number = 8;
bm(8).t1 = '2014-feb-15:10:24:42';
bm(8).t2 = '2014-feb-15:18:06:27';
bm(8).dk = 58;
bm(8).sched = '7_ffflat_dslmast_31h_006.sch';
bm(8).filename = '20140215_102917';
bm(8).notes = '';

bm(9).number = 9;
bm(9).t1 = '2014-feb-15:18:06:28';
bm(9).t2 = '2014-feb-16:01:46:56';
bm(9).dk = 130;
bm(9).sched = '7_ffflat_dslmast_31i_006.sch';
bm(9).filename = '20140215_181103';
bm(9).notes = '';

bm(10).number = 10;
bm(10).t1 = '2014-feb-16:01:46:57';
bm(10).t2 = '2014-feb-16:09:27:26';
bm(10).dk = 202;
bm(10).sched = '7_ffflat_dslmast_31j_006.sch';
bm(10).filename = '20140216_015133';
bm(10).notes = '';

bm(11).number = 11;
bm(11).t1 = '2014-feb-16:15:32:43';
bm(11).t2 = '2014-feb-16:23:14:11';
bm(11).dk = -50;
bm(11).sched = '7_ffflat_dslmast_31b_007.sch';
bm(11).filename = '20140216_153718';
bm(11).notes = 'Were the squids tuned?';

bm(12).number = 12;
bm(12).t1 = '2014-feb-17:00:14:54';
bm(12).t2 = '2014-feb-17:07:54:32';
bm(12).dk = -122;
bm(12).sched = '7_ffflat_dslmast_31a_007.sch';
bm(12).filename = '20140217_001454';
bm(12).notes = 'Were the squids tuned?';

bm(13).number = 13;
bm(13).t1 = '2014-feb-17:07:54:32';
bm(13).t2 = '2014-feb-17:15:35:47';
bm(13).dk = 22;
bm(13).sched = '7_ffflat_dslmast_31c_007.sch';
bm(13).filename = '20140217_075908';
bm(13).notes = '';

bm(14).number = 14;
bm(14).t1 = '2014-feb-17:15:35:47';
bm(14).t2 = '2014-feb-17:23:17:32';
bm(14).dk = 94;
bm(14).sched = '7_ffflat_dslmast_31d_007.sch';
bm(14).filename = '20140217_153345';
bm(14).notes = '';

bm(15).number = 15;
bm(15).t1 = '2014-feb-18:08:08:50';
bm(15).t2 = '2014-feb-18:15:48:28';
bm(15).dk = 166;
bm(15).sched = '7_ffflat_dslmast_31e_007.sch';
bm(15).filename = '20140218_081325';
bm(15).notes = '';

bm(16).number = 16;
bm(16).t1 = '2014-feb-18:15:48:28';
bm(16).t2 = '2014-feb-18:23:29:21';
bm(16).dk = -86;
bm(16).sched = '7_ffflat_dslmast_31f_007.sch';
bm(16).filename = '20140218_154747';
bm(16).notes = '';

bm(17).number = 17;
bm(17).t1 = '2014-feb-18:23:29:22';
bm(17).t2 = '2014-feb-19:07:09:50';
bm(17).dk = -14;
bm(17).sched = '7_ffflat_dslmast_31g_007.sch';
bm(17).filename = '20140218_233357';
bm(17).notes = '';

bm(18).number = 18;
bm(18).t1 = '2014-feb-19:07:09:51';
bm(18).t2 = '2014-feb-19:14:51:36';
bm(18).dk = 58;
bm(18).sched = '7_ffflat_dslmast_31h_007.sch';
bm(18).filename = '20140219_071426';
bm(18).notes = '';

bm(19).number = 19;
bm(19).t1 = '2014-feb-19:20:39:17';
bm(19).t2 = '2014-feb-20:04:19:14';
bm(19).dk = 130;
bm(19).sched = '7_ffflat_dslmast_31i_007.sch';
bm(19).filename = '20140219_204347';
bm(19).notes = '';

bm(20).number = 20;
bm(20).t1 = '2014-feb-20:04:19:15';
bm(20).t2 = '2014-feb-20:11:59:39';
bm(20).dk = 202;
bm(20).sched = '7_ffflat_dslmast_31j_007.sch';
bm(20).filename = '20140220_042350';
bm(20).notes = '';

bm(21).number = 21;
bm(21).t1 = '2014-feb-20:11:59:40';
bm(21).t2 = '2014-feb-20:19:40:34';
bm(21).dk = -50;
bm(21).sched = '7_ffflat_dslmast_31b_008.sch';
bm(21).filename = '20140220_120415';
bm(21).notes = '';

bm(22).number = 22;
bm(22).t1 = '2014-feb-20:19:40:35';
bm(22).t2 = '2014-feb-21:03:21:23';
bm(22).dk = -122;
bm(22).sched = '7_ffflat_dslmast_31a_008.sch';
bm(22).filename = '20140220_194510';
bm(22).notes = '';

bm(23).number = 23;
bm(23).t1 = '2014-feb-21:09:19:33';
bm(23).t2 = '2014-feb-21:17:00:08';
bm(23).dk = 22;
bm(23).sched = '7_ffflat_dslmast_31c_008.sch';
bm(23).filename = '20140221_092408';
bm(23).notes = '';

bm(24).number = 24;
bm(24).t1 = '2014-feb-21:17:00:09';
bm(24).t2 = '2014-feb-22:00:41:54';
bm(24).dk = 94;
bm(24).sched = '7_ffflat_dslmast_31d_008.sch';
bm(24).filename = '20140221_170444';
bm(24).notes = '';

bm(25).number = 25;
bm(25).t1 = '2014-feb-22:02:44:56';
bm(25).t2 = '2014-feb-22:10:26:08';
bm(25).dk = -50;
bm(25).sched = '7_ffflat_dslmast_31b_009.sch';
bm(25).filename = '20140222_024932';
bm(25).notes = '';

bm(26).number = 26;
bm(26).t1 = '2014-feb-22:10:26:09';
bm(26).t2 = '2014-feb-22:18:06:46';
bm(26).dk = -122;
bm(26).sched = '7_ffflat_dslmast_31a_009.sch';
bm(26).filename = '20140222_103044';
bm(26).notes = '';

bm(27).number = 27;
bm(27).t1 = '2014-feb-22:23:50:45';
bm(27).t2 = '2014-feb-23:07:31:40';
bm(27).dk = 22;
bm(27).sched = '7_ffflat_dslmast_31c_009.sch';
bm(27).filename = '20140222_235514';
bm(27).notes = '';

bm(28).number = 28;
bm(28).t1 = '2014-feb-23:07:31:41';
bm(28).t2 = '2014-feb-23:15:13:25';
bm(28).dk = 94;
bm(28).sched = '7_ffflat_dslmast_31d_009.sch';
bm(28).filename = '20140223_073616';
bm(28).notes = '';

bm(29).number = 29;
bm(29).t1 = '2014-feb-23:15:13:26';
bm(29).t2 = '2014-feb-23:22:54:05';
bm(29).dk = 166;
bm(29).sched = '7_ffflat_dslmast_31e_009.sch';
bm(29).filename = '20140223_151801';
bm(29).notes = '';

bm(30).number = 30;
bm(30).t1 = '2014-feb-23:22:54:06';
bm(30).t2 = '2014-feb-24:06:35:15';
bm(30).dk = -86;
bm(30).sched = '7_ffflat_dslmast_31f_009.sch';
bm(30).filename = '20140223_225841';
bm(30).notes = '';

bm(31).number = 31;
bm(31).t1 = '2014-feb-24:11:55:51';
bm(31).t2 = '2014-feb-24:19:35:46';
bm(31).dk = -14;
bm(31).sched = '7_ffflat_dslmast_31g_009.sch';
bm(31).filename = '20140224_120019';
bm(31).notes = '';

bm(32).number = 32;
bm(32).t1 = '2014-feb-24:19:35:47';
bm(32).t2 = '2014-feb-25:03:17:31';
bm(32).dk = 58;
bm(32).sched = '7_ffflat_dslmast_31h_009.sch';
bm(32).filename = '20140224_194022';
bm(32).notes = '';

bm(33).number = 33;
bm(33).t1 = '2014-feb-25:03:17:32';
bm(33).t2 = '2014-feb-25:10:57:59';
bm(33).dk = 130;
bm(33).sched = '7_ffflat_dslmast_31i_009.sch';
bm(33).filename = '20140225_032207';
bm(33).notes = '';

bm(34).number = 34;
bm(34).t1 = '2014-feb-25:10:57:59';
bm(34).t2 = '2014-feb-25:18:38:43';
bm(34).dk = 202;
bm(34).sched = '7_ffflat_dslmast_31j_009.sch';
bm(34).filename = '20140225_110234';
bm(34).notes = '';

bm(35).number = 35;
bm(35).t1 = '2014-feb-26:03:08:51';
bm(35).t2 = '2014-feb-26:10:48:49';
bm(35).dk = -122;
bm(35).sched = '7_ffflat_dslmast_31a_010.sch';
bm(35).filename = '20140226_031317';
bm(35).notes = '';

bm(36).number = 36;
bm(36).t1 = '2014-feb-26:10:48:50';
bm(36).t2 = '2014-feb-26:18:29:08';
bm(36).dk = -50;
bm(36).sched = '7_ffflat_dslmast_31b_010.sch';
bm(36).filename = '20140226_105325';
bm(36).notes = '';

bm(37).number = 37;
bm(37).t1 = '2014-feb-26:18:29:09';
bm(37).t2 = '2014-feb-27:02:09:58';
bm(37).dk = 22;
bm(37).sched = '7_ffflat_dslmast_31c_010.sch';
bm(37).filename = '20140226_183344';
bm(37).notes = '';

bm(38).number = 38;
bm(38).t1 = '2014-feb-27:02:09:59';
bm(38).t2 = '2014-feb-27:09:51:45';
bm(38).dk = 94;
bm(38).sched = '7_ffflat_dslmast_31d_010.sch';
bm(38).filename = '20140227_021434';
bm(38).notes = '';

bm(39).number = 39;
bm(39).t1 = '2014-feb-27:19:31:37';
bm(39).t2 = '2014-feb-28:03:11:33';
bm(39).dk = 166;
bm(39).sched = '7_ffflat_dslmast_31e_010.sch';
bm(39).filename = '20140227_193606';
bm(39).notes = '';

bm(40).number = 40;
bm(40).t1 = '2014-feb-28:03:11:34';
bm(40).t2 = '2014-feb-28:10:52:31';
bm(40).dk = -86;
bm(40).sched = '7_ffflat_dslmast_31f_010.sch';
bm(40).filename = '20140228_031609';
bm(40).notes = '';

bm(41).number = 41;
bm(41).t1 = '2014-feb-28:10:52:32';
bm(41).t2 = '2014-feb-28:18:32:35';
bm(41).dk = -14;
bm(41).sched = '7_ffflat_dslmast_31g_010.sch';
bm(41).filename = '20140228_105707';
bm(41).notes = '';

bm(42).number = 42;
bm(42).t1 = '2014-feb-28:18:32:35';
bm(42).t2 = '2014-mar-01:02:14:20';
bm(42).dk = 58;
bm(42).sched = '7_ffflat_dslmast_31h_010.sch';
bm(42).filename = '20140228_183144';
bm(42).notes = '';

bm(43).number = 43;
bm(43).t1 = '2014-mar-01:07:35:37';
bm(43).t2 = '2014-mar-01:15:15:35';
bm(43).dk = 130;
bm(43).sched = '7_ffflat_dslmast_31i_010.sch';
bm(43).filename = '20140301_074009';
bm(43).notes = '';

bm(44).number = 44;
bm(44).t1 = '2014-mar-01:15:15:36';
bm(44).t2 = '2014-mar-01:22:55:38';
bm(44).dk = 202;
bm(44).sched = '7_ffflat_dslmast_31j_010.sch';
bm(44).filename = '20140301_152011';
bm(44).notes = '';

bm(45).number = 45;
bm(45).t1 = '2014-mar-01:22:55:39';
bm(45).t2 = '2014-mar-02:06:36:44';
bm(45).dk = 58;
bm(45).sched = '7_ffflat_dslmast_31h_011.sch';
bm(45).filename = '20140301_230014';
bm(45).notes = '';

bm(46).number = 46;
bm(46).t1 = '2014-mar-02:06:36:45';
bm(46).t2 = '2014-mar-02:14:19:28';
bm(46).dk = -14;
bm(46).sched = '7_ffflat_dslmast_31g_011.sch';
bm(46).filename = '20140302_064502';
bm(46).notes = '';

bm(47).number = 47;
bm(47).t1 = '2014-mar-02:23:10:10';
bm(47).t2 = '2014-mar-03:06:51:29';
bm(47).dk = 202;
bm(47).sched = '7_ffflat_dslmast_31j_011.sch';
bm(47).filename = '20140302_231445';
bm(47).notes = '';

bm(48).number = 48;
bm(48).t1 = '2014-mar-04:01:59:41';
bm(48).t2 = '2014-mar-04:09:39:21';
bm(48).dk = -122;
bm(48).sched = '7_ffflat_dslmast_31a_013.sch';
bm(48).filename = '20140304_020411';
bm(48).notes = '';

bm(49).number = 49;
bm(49).t1 = '2014-mar-04:09:39:22';
bm(49).t2 = '2014-mar-04:17:19:48';
bm(49).dk = -50;
bm(49).sched = '7_ffflat_dslmast_31b_013.sch';
bm(49).filename = '20140304_094358';
bm(49).notes = '';

bm(50).number = 50;
bm(50).t1 = '2014-mar-04:17:19:49';
bm(50).t2 = '2014-mar-05:01:00:39';
bm(50).dk = 22;
bm(50).sched = '7_ffflat_dslmast_31c_013.sch';
bm(50).filename = '20140304_172424';
bm(50).notes = '';

bm(51).number = 51;
bm(51).t1 = '2014-mar-05:01:00:40';
bm(51).t2 = '2014-mar-05:08:42:29';
bm(51).dk = 94;
bm(51).sched = '7_ffflat_dslmast_31d_013.sch';
bm(51).filename = '20140305_010515';
bm(51).notes = '';

bm(52).number = 52;
bm(52).t1 = '2014-mar-05:14:15:22';
bm(52).t2 = '2014-mar-05:21:57:12';
bm(52).dk = 58;
bm(52).sched = '7_ffflat_dslmast_31h_013.sch';
bm(52).filename = '20140305_141953';
bm(52).notes = '';

bm(53).number = 53;
bm(53).t1 = '2014-mar-05:21:57:13';
bm(53).t2 = '2014-mar-06:05:38:59';
bm(53).dk = -14;
bm(53).sched = '7_ffflat_dslmast_31g_013.sch';
bm(53).filename = '20140305_220148';
bm(53).notes = '';

bm(54).number = 54;
bm(54).t1 = '2014-mar-06:05:39:00';
bm(54).t2 = '2014-mar-06:13:19:54';
bm(54).dk = 166;
bm(54).sched = '7_ffflat_dslmast_31e_013.sch';
bm(54).filename = '20140306_054335';
bm(54).notes = '';

bm(55).number = 55;
bm(55).t1 = '2014-mar-06:13:19:55';
bm(55).t2 = '2014-mar-06:21:00:14';
bm(55).dk = -86;
bm(55).sched = '7_ffflat_dslmast_31f_013.sch';
bm(55).filename = '20140306_132430';
bm(55).notes = '';

bm(56).number = 56;
bm(56).t1 = '2014-mar-07:02:12:42';
bm(56).t2 = '2014-mar-07:09:53:14';
bm(56).dk = 130;
bm(56).sched = '7_ffflat_dslmast_31i_013.sch';
bm(56).filename = '20140307_021713';
bm(56).notes = '';

bm(57).number = 57;
bm(57).t1 = '2014-mar-07:09:53:15';
bm(57).t2 = '2014-mar-07:17:33:40';
bm(57).dk = 202;
bm(57).sched = '7_ffflat_dslmast_31j_013.sch';
bm(57).filename = '20140307_095750';
bm(57).notes = '';


return


