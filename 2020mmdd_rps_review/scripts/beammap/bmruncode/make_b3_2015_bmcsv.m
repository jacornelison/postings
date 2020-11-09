function make_b3_2015_bmccsv()
% Script to generate a nice archival .csv for the 2015 BICEP3 beam maps
% Script is now obsolete
% Please use make_bmcsv which is more generic and yields the same result with 10x less work

disp('Script is now obsolete.')
disp(' Please use make_bmcsv with appropriate exp and year. It is more general  and yields the same result with 10x less work')
return 

% Don't overwrite a different experiment!
if ~strcmp(get_experiment_name(),'bicep3');
  error('Get into the right folder for BICEP3!')
end

savefile = 'beammaps/bmrunlist_2015.csv';

bm = bm_info();

for i = 1:length(bm)
  p.number(i) = bm(i).number;
  p.t1{i} = bm(i).t1;
  p.t2{i} = bm(i).t2;
  p.dk(i) = bm(i).dk;
  p.sched{i} = bm(i).sched;
  p.filename{i} = bm(i).filename;
  p.notes{i} = bm(i).notes;
  
end

k.comments{1} = ['# BEAM MAPPING RUNS: BICEP3 2015'];
k.fields = {'number','t1','t2','dk','sched','filename','notes'};
k.formats = ...
    {'integer','string','string','double','string','string','string'};
k.units = ...
    {'unitless','time','time','deg','unitless','time','unitless'};

ParameterWrite(savefile,p,k);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bm = bm_info()

bm(1).number = 1;
bm(1).t1 = '2015-feb-10:03:06:29';
bm(1).t2 = '2015-feb-10:03:54:10';
bm(1).dk = 0;
bm(1).sched = '7_ffflat_mapo_12_001.sch';
bm(1).filename = '';
bm(1).notes = '';

bm(2).number = 2;
bm(2).t1 = '2015-feb-10:04:34:17';
bm(2).t2 = '2015-feb-10:06:44:23';
bm(2).dk = 0;
bm(2).sched = '7_ffflat_mapo_13_001.sch';
bm(2).filename = '';
bm(2).notes = '';

bm(3).number = 3;
bm(3).t1 = '2015-feb-11:23:06:55';
bm(3).t2 = '2015-feb-12:06:19:13';
bm(3).dk = 0;
bm(3).sched = '7_ffflat_mapo_33_dk090_001.sch';
bm(3).filename = '';
bm(3).notes = '';

bm(4).number = 4;
bm(4).t1 = '2015-feb-12:06:53:12';
bm(4).t2 = '2015-feb-12:06:42:00';
bm(4).dk = 90;
bm(4).sched = '7_ffflat_mapo_33_dk090_002.sch';
bm(4).filename = '';
bm(4).notes = '';

bm(5).number = 5;
bm(5).t1 = '2015-feb-13:01:12:32';
bm(5).t2 = '2015-feb-13:09:54:00';
bm(5).dk = 0;
bm(5).sched = '7_ffflat_mapo_34_dk000_001.sch';
bm(5).filename = '';
bm(5).notes = '';

bm(6).number = 6;
bm(6).t1 = '2015-feb-13:10:11:12';
bm(6).t2 = '2015-feb-13:19:22:15';
bm(6).dk = 90;
bm(6).sched = '7_ffflat_mapo_34_dk090_001.sch';
bm(6).filename = '';
bm(6).notes = '';

bm(7).number = 7;
bm(7).t1 = '2015-feb-14:06:45:46';
bm(7).t2 = '2015-feb-14:14:29:28';
bm(7).dk = 90;
bm(7).sched = '7_ffflat_mapo_35_dk090_001.sch';
bm(7).filename = '20150214_064602';
bm(7).notes = '';

bm(8).number = 8;
bm(8).t1 = '2015-feb-15:06:51:13';
bm(8).t2 = '2015-feb-15:14:34:40';
bm(8).dk = 45;
bm(8).sched = '7_ffflat_mapo_35_dk045_001.sch';
bm(8).filename = '20150215_065128';
bm(8).notes = '';

bm(9).number = 9;
bm(9).t1 = '2015-feb-16:05:32:14';
bm(9).t2 = '2015-feb-16:13:15:52';
bm(9).dk = 315;
bm(9).sched = '7_ffflat_mapo_35_dk315_001.sch';
bm(9).filename = '20150216_053229';
bm(9).notes = '';

bm(10).number = 10;
bm(10).t1 = '2015-feb-17:05:51:20';
bm(10).t2 = '2015-feb-17:12:52:47';
bm(10).dk = -45;
bm(10).sched = '7_ffflat_mapo_35_dk270_001.sch';
bm(10).filename = '20150217_055136';
bm(10).notes = 'Check dk';

bm(11).number = 11;
bm(11).t1 = '2015-feb-18:05:33:10';
bm(11).t2 = '2015-feb-18:13:17:01';
bm(11).dk = -90;
bm(11).sched = '7_ffflat_mapo_35_dk270_002.sch';
bm(11).filename = '20150218_053325';
bm(11).notes = '';

bm(12).number = 12;
bm(12).t1 = '2015-feb-19:05:33:48';
bm(12).t2 = '2015-feb-19:13:17:17';
bm(12).dk = 135;
bm(12).sched = '7_ffflat_mapo_35_dk135_001.sch';
bm(12).filename = '';
bm(12).notes = '';

bm(13).number = 13;
bm(13).t1 = '2015-feb-20:05:47:16';
bm(13).t2 = '2015-feb-20:13:30:53';
bm(13).dk = 135;
bm(13).sched = '7_ffflat_mapo_35_dk135_002.sch';
bm(13).filename = '';
bm(13).notes = '';

bm(14).number = 14;
bm(14).t1 = '2015-feb-21:05:35:55';
bm(14).t2 = '2015-feb-21:13:18:15';
bm(14).dk = 157.5;
bm(14).sched = '7_ffflat_mapo_35_dk157-5_001.sch';
bm(14).filename = '';
bm(14).notes = '';

bm(15).number = 15;
bm(15).t1 = '2015-feb-22:05:34:34';
bm(15).t2 = '2015-feb-22:13:18:15';
bm(15).dk = 112.5;
bm(15).sched = '7_ffflat_mapo_35_dk112p5_001.sch';
bm(15).filename = '';
bm(15).notes = '';

bm(16).number = 16;
bm(16).t1 = '2015-feb-23:05:38:18';
bm(16).t2 = '2015-feb-23:13:21:46';
bm(16).dk = 67.5;
bm(16).sched = '7_ffflat_mapo_35_dk067p5_001.sch';
bm(16).filename = '';
bm(16).notes = '';

bm(17).number = 17;
bm(17).t1 = '2015-feb-26:06:05:09';
bm(17).t2 = '2015-feb-26:13:48:43';
bm(17).dk = 22.5;
bm(17).sched = '7_ffflat_mapo_35_dk022p5_xxx.sch';
bm(17).filename = '20150226_060524';
bm(17).notes = '';

bm(18).number = 18;
bm(18).t1 = '2015-feb-27:05:39:04';
bm(18).t2 = '2015-feb-27:13:22:38';
bm(18).dk = 337.5;
bm(18).sched = '7_ffflat_mapo_35_dk337p5_001.sch';
bm(18).filename = '20150227_053919';
bm(18).notes = '';

bm(19).number = 19;
bm(19).t1 = '2015-mar-01:10:05:56';
bm(19).t2 = '2015-mar-01:17:49:39';
bm(19).dk = 292.5;
bm(19).sched = '7_ffflat_mapo_35_dk292p5_001.sch';
bm(19).filename = '';
bm(19).notes = '';

bm(20).number = 20;
bm(20).t1 = '2015-mar-02:05:32:05';
bm(20).t2 = '2015-mar-02:13:15:39';
bm(20).dk = 157.5;
bm(20).sched = '7_ffflat_mapo_35_dk157p5_002.sch';
bm(20).filename = '';
bm(20).notes = '';

bm(21).number = 21;
bm(21).t1 = '2015-mar-03:05:33:17';
bm(21).t2 = '2015-mar-03:13:16:59';
bm(21).dk = 135;
bm(21).sched = '7_ffflat_mapo_35_dk135_003.sch';
bm(21).filename = '';
bm(21).notes = '';

bm(22).number = 22;
bm(22).t1 = '2015-mar-04:05:35:19';
bm(22).t2 = '2015-mar-04:13:18:48';
bm(22).dk = 112.5;
bm(22).sched = '7_ffflat_mapo_35_dk112p5_002.sch';
bm(22).filename = '';
bm(22).notes = '';

bm(23).number = 23;
bm(23).t1 = '2015-mar-05:05:37:42';
bm(23).t2 = '2015-mar-05:13:21:29';
bm(23).dk = 90;
bm(23).sched = '7_ffflat_mapo_35_dk090_002.sch';
bm(23).filename = '';
bm(23).notes = '';

bm(24).number = 24;
bm(24).t1 = '2015-mar-06:05:32:40';
bm(24).t2 = '2015-mar-06:13:16:08';
bm(24).dk = 67.5;
bm(24).sched = '7_ffflat_mapo_35_dk067p5_002.sch';
bm(24).filename = '';
bm(24).notes = '';

bm(25).number = 25;
bm(25).t1 = '2015-mar-06:20:03:39';
bm(25).t2 = '2015-mar-07:00:43:46';
bm(25).dk = 0;
bm(25).sched = '7_ffflat_mapo_37b_dk000_001.sch';
bm(25).filename = '20150306_200354';
bm(25).notes = '';

bm(26).number = 26;
bm(26).t1 = '2015-mar-07:08:01:52';
bm(26).t2 = '2015-mar-07:12:42:00';
bm(26).dk = 0;
bm(26).sched = '7_ffflat_mapo_37c_dk000_001.sch';
bm(26).filename = '20150307_080208';
bm(26).notes = '';

bm(27).number = 27;
bm(27).t1 = '2015-mar-07:19:18:32';
bm(27).t2 = '2015-mar-07:23:58:40';
bm(27).dk = 90;
bm(27).sched = '7_ffflat_mapo_37b_dk090_001.sch';
bm(27).filename = '20150307_191848';
bm(27).notes = '';

bm(28).number = 28;
bm(28).t1 = '2015-mar-08:06:33:38';
bm(28).t2 = '2015-mar-08:11:13:45';
bm(28).dk = 90;
bm(28).sched = '7_ffflat_mapo_37c_dk090_001.sch';
bm(28).filename = '';
bm(28).notes = '';

bm(29).number = 29;
bm(29).t1 = '2015-mar-10:05:27:30';
bm(29).t2 = '2015-mar-10:10:07:38';
bm(29).dk = 0;
bm(29).sched = '2_ffflat_mapo_37a_dk000_001.sch';
bm(29).filename = '20150310_052745';
bm(29).notes = '';

bm(30).number = 30;
bm(30).t1 = '2015-mar-10:16:43:00';
bm(30).t2 = '2015-mar-10:21:23:19';
bm(30).dk = 0;
bm(30).sched = '4_ffflat_mapo_37d_dk000_001.sch';
bm(30).filename = '20150310_164315';
bm(30).notes = '';

bm(31).number = 31;
bm(31).t1 = '2015-mar-11:05:31:17';
bm(31).t2 = '2015-mar-11:10:11:25';
bm(31).dk = 90;
bm(31).sched = '2_ffflat_mapo_37a_dk090_001.sch';
bm(31).filename = '';
bm(31).notes = '';

bm(32).number = 32;
bm(32).t1 = '2015-mar-11:16:45:48';
bm(32).t2 = '2015-mar-11:21:25:55';
bm(32).dk = 90;
bm(32).sched = '4_ffflat_mapo_37d_dk090_001.sch';
bm(32).filename = '';
bm(32).notes = '';

bm(33).number = 33;
bm(33).t1 = '2015-mar-12:05:36:50';
bm(33).t2 = '2015-mar-12:10:17:05';
bm(33).dk = 45;
bm(33).sched = '2_ffflat_mapo_37a_dk045_001.sch';
bm(33).filename = '';
bm(33).notes = '';

bm(34).number = 34;
bm(34).t1 = '2015-mar-12:17:02:18';
bm(34).t2 = '2015-mar-12:21:42:26';
bm(34).dk = 45;
bm(34).sched = '4_ffflat_mapo_37b_dk045_001.sch';
bm(34).filename = '';
bm(34).notes = '';

bm(35).number = 35;
bm(35).t1 = '2015-mar-13:05:31:08';
bm(35).t2 = '2015-mar-13:10:11:16';
bm(35).dk = 45;
bm(35).sched = '2_ffflat_mapo_37c_dk045_001.sch';
bm(35).filename = '';
bm(35).notes = '';

bm(36).number = 36;
bm(36).t1 = '2015-mar-13:16:46:19';
bm(36).t2 = '2015-mar-13:21:26:27';
bm(36).dk = 45;
bm(36).sched = '4_ffflat_mapo_37d_dk045_001.sch';
bm(36).filename = '';
bm(36).notes = '';

bm(37).number = 37;
bm(37).t1 = '2015-mar-14:10:36:00';
bm(37).t2 = '2015-mar-14:15:16:17';
bm(37).dk = 135;
bm(37).sched = '2_ffflat_mapo_37a_dk135_002.sch';
bm(37).filename = '';
bm(37).notes = '';

bm(38).number = 38;
bm(38).t1 = '2015-mar-14:21:48:27';
bm(38).t2 = '2015-mar-15:02:28:43';
bm(38).dk = 135;
bm(38).sched = '4_ffflat_mapo_37b_dk135_002.sch';
bm(38).filename = '';
bm(38).notes = '';

bm(39).number = 39;
bm(39).t1 = '2015-mar-15:09:03:59';
bm(39).t2 = '2015-mar-15:13:44:09';
bm(39).dk = 135;
bm(39).sched = '4_ffflat_mapo_38c_dk135_001.sch';
bm(39).filename = '';
bm(39).notes = '';

bm(40).number = 40;
bm(40).t1 = '2015-mar-15:20;24:10';
bm(40).t2 = '2015-mar-16:01:04:20';
bm(40).dk = 135;
bm(40).sched = '4_ffflat_mapo_38d_dk135_001.sch';
bm(40).filename = '';
bm(40).notes = '';

bm(41).number = 41;
bm(41).t1 = '2015-mar-16:07:38:32';
bm(41).t2 = '2015-mar-16:12:18:41';
bm(41).dk = 0;
bm(41).sched = '6_ffflat_mapo_38a_dk000_001.sch';
bm(41).filename = '';
bm(41).notes = '';

bm(42).number = 42;
bm(42).t1 = '2015-mar-16:18:51:00';
bm(42).t2 = '2015-mar-16:23:31:14';
bm(42).dk = 0;
bm(42).sched = '8_ffflat_mapo_38b_dk000_001.sch';
bm(42).filename = '';
bm(42).notes = '';

bm(43).number = 43;
bm(43).t1 = '2015-mar-17:12:34:56';
bm(43).t2 = '2015-mar-17:17:15:04';
bm(43).dk = 0;
bm(43).sched = '2_ffflat_mapo_38a_dk000_002.sch';
bm(43).filename = '20150317_123511';
bm(43).notes = '';

bm(44).number = 44;
bm(44).t1 = '2015-mar-18:00:02:01';
bm(44).t2 = '2015-mar-18:04:42:08';
bm(44).dk = 0;
bm(44).sched = '4_ffflat_mapo_38b_dk000_002.sch';
bm(44).filename = '20150317_123511';
bm(44).notes = '';

bm(45).number = 45;
bm(45).t1 = '2015-mar-18:11:12:41';
bm(45).t2 = '2015-mar-18:15:52:49';
bm(45).dk = 0;
bm(45).sched = '6_ffflat_mapo_38c_dk000_001.sch';
bm(45).filename = '20150318_111257';
bm(45).notes = '';

bm(46).number = 46;
bm(46).t1 = '2015-mar-18:22:23:27';
bm(46).t2 = '2015-mar-19:03:03:34';
bm(46).dk = 0;
bm(46).sched = '8_ffflat_mapo_38d_dk000_001.sch';
bm(46).filename = '20150318_222342';
bm(46).notes = '';

bm(47).number = 47;
bm(47).t1 = '2015-mar-19:09:40:31';
bm(47).t2 = '2015-mar-19:14:20:39';
bm(47).dk = 45;
bm(47).sched = '2_ffflat_mapo_38a_dk045_001.sch';
bm(47).filename = '20150319_094047';
bm(47).notes = '';

bm(48).number = 48;
bm(48).t1 = '2015-mar-19:20:53:58';
bm(48).t2 = '2015-mar-20:01:34:06';
bm(48).dk = 45;
bm(48).sched = '4_ffflat_mapo_38b_dk045_001.sch';
bm(48).filename = '20150319_205414';
bm(48).notes = '';

bm(49).number = 49;
bm(49).t1 = '2015-mar-20:08:06:16';
bm(49).t2 = '2015-mar-20:12:46:24';
bm(49).dk = 45;
bm(49).sched = '6_ffflat_mapo_38c_dk045_001.sch';
bm(49).filename = '20150320_080631';
bm(49).notes = '';

bm(50).number = 50;
bm(50).t1 = '2015-mar-20:19:23:32';
bm(50).t2 = '2015-mar-21:00:03:53';
bm(50).dk = 45;
bm(50).sched = '8_ffflat_mapo_38d_dk045_001.sch';
bm(50).filename = '20150320_192348';
bm(50).notes = 'mce2 offline between 23:40:36 and 23:41:37';

bm(51).number = 51;
bm(51).t1 = '2015-mar-21:07:01:24';
bm(51).t2 = '2015-mar-21:11:41:38';
bm(51).dk = 90;
bm(51).sched = '2_ffflat_mapo_38a_dk090_001.sch';
bm(51).filename = '20150321_070139';
bm(51).notes = '';

bm(52).number = 52;
bm(52).t1 = '2015-mar-21:18:35:12';
bm(52).t2 = '2015-mar-21:23:15:19';
bm(52).dk = 90;
bm(52).sched = '4_ffflat_mapo_38b_dk090_001.sch';
bm(52).filename = '20150321_183527';
bm(52).notes = '';

bm(53).number = 53;
bm(53).t1 = '2015-mar-22:06:07:30';
bm(53).t2 = '2015-mar-22:10:47:38';
bm(53).dk = 90;
bm(53).sched = '6_ffflat_mapo_38c_dk090_001.sch';
bm(53).filename = '20150322_060746';
bm(53).notes = '';

bm(54).number = 54;
bm(54).t1 = '2015-mar-23:02:16:21';
bm(54).t2 = '2015-mar-23:06:58:10';
bm(54).dk = 45;
bm(54).sched = '2_ffflat_mapo_38d_dk045_002.sch';
bm(54).filename = '20150323_021636';
bm(54).notes = '';

bm(55).number = 55;
bm(55).t1 = '2015-mar-24:00:54:54';
bm(55).t2 = '2015-mar-24:05:36:42';
bm(55).dk = 90;
bm(55).sched = '2_ffflat_mapo_38a_dk090_002.sch';
bm(55).filename = '';
bm(55).notes = '';

bm(56).number = 56;
bm(56).t1 = '2015-mar-24:12:07:16';
bm(56).t2 = '2015-mar-24:06:49:05';
bm(56).dk = 90;
bm(56).sched = '4_ffflat_mapo_38b_dk090_002.sch';
bm(56).filename = '';
bm(56).notes = '';

bm(57).number = 57;
bm(57).t1 = '2015-mar-24:23:21:41';
bm(57).t2 = '2015-mar-25:04:03:34';
bm(57).dk = 90;
bm(57).sched = '6_ffflat_mapo_38c_dk090_002.sch';
bm(57).filename = '';
bm(57).notes = '';

bm(58).number = 58;
bm(58).t1 = '2015-mar-25:10:38:58';
bm(58).t2 = '2015-mar-25:15:20:46';
bm(58).dk = 90;
bm(58).sched = '8_ffflat_mapo_38d_dk090_001.sch';
bm(58).filename = '';
bm(58).notes = 'Chop ref lost in schedule';

bm(59).number = 59;
bm(59).t1 = '2015-mar-25:21:53:37';
bm(59).t2 = '2015-mar-26:04:26:51';
bm(59).dk = 135;
bm(59).sched = '2_ffflat_mapo_39a_dk135_001.sch';
bm(59).filename = '';
bm(59).notes = 'No chop ref';

bm(60).number = 60;
bm(60).t1 = '2015-mar-26:11:03:51';
bm(60).t2 = '2015-mar-26:18:04:55';
bm(60).dk = 0;
bm(60).sched = '2_ffflat_mapo_40b_dk000_001.sch';
bm(60).filename = '';
bm(60).notes = 'No chop ref';

bm(61).number = 61;
bm(61).t1 = '2015-mar-27:01:05:56';
bm(61).t2 = '2015-mar-27:04:54:47';
bm(61).dk = 0;
bm(61).sched = '4_ffflat_mapo_40a_dk000_001.sch';
bm(61).filename = '';
bm(61).notes = 'No chop ref';

bm(62).number = 62;
bm(62).t1 = '2015-mar-27:11:56:31';
bm(62).t2 = '2015-mar-27:18:29:38';
bm(62).dk = 0;
bm(62).sched = '4_ffflat_mapo_40a_dk000_002.sch';
bm(62).filename = '';
bm(62).notes = 'Chop ref back? (check)';

bm(63).number = 63;
bm(63).t1 = '2015-mar-28:01:01:09';
bm(63).t2 = '2015-mar-28:07:34:14';
bm(63).dk = 45;
bm(63).sched = '6_ffflat_mapo_40a_dk045_001.sch';
bm(63).filename = '20150328_010124';
bm(63).notes = '';

bm(64).number = 64;
bm(64).t1 = '2015-mar-28:14:25:13';
bm(64).t2 = '2015-mar-28:20:58:41';
bm(64).dk = 45;
bm(64).sched = '8_ffflat_mapo_40b_dk045_001.sch';
bm(64).filename = '20150328_142529';
bm(64).notes = '';

bm(65).number = 65;
bm(65).t1 = '2015-mar-29:03:29:55';
bm(65).t2 = '2015-mar-29:10:03:03';
bm(65).dk = 90;
bm(65).sched = '2_ffflat_mapo_40a_dk090_001.sch';
bm(65).filename = '20150329_033010';
bm(65).notes = '';

bm(66).number = 66;
bm(66).t1 = '2015-mar-29:16:34:41';
bm(66).t2 = '2015-mar-29:23:07:49';
bm(66).dk = 90;
bm(66).sched = '4_ffflat_mapo_40b_dk090_001.sch';
bm(66).filename = '20150329_163456';
bm(66).notes = '';

bm(67).number = 67;
bm(67).t1 = '2015-mar-30:05:39:26';
bm(67).t2 = '2015-mar-30:09:22:43';
bm(67).dk = 135;
bm(67).sched = '6_ffflat_mapo_40a_dk135_001.sch';
bm(67).filename = '';
bm(67).notes = 'mce0 lost between 08:25:53 and 08:26:54';

bm(68).number = 68;
bm(68).t1 = '2015-mar-30:16:07:21';
bm(68).t2 = '2015-mar-30:22:40:32';
bm(68).dk = 135;
bm(68).sched = '6_ffflat_mapo_40a_dk135_002.sch';
bm(68).filename = '20150330_160736';
bm(68).notes = '';

bm(69).number = 69;
bm(69).t1 = '2015-mar-31:05:11:55';
bm(69).t2 = '2015-mar-31:11:45:37';
bm(69).dk = 135;
bm(69).sched = '8_ffflat_mapo_40b_dk135_001.sch';
bm(69).filename = '20150331_051211';
bm(69).notes = '';

bm(70).number = 70;
bm(70).t1 = '2015-mar-31:18:26:39';
bm(70).t2 = '2015-apr-01:00:59:49';
bm(70).dk = 0;
bm(70).sched = '2_ffflat_mapo_40b_dk000_002.sch';
bm(70).filename = '20150331_182654';
bm(70).notes = '';

bm(71).number = 71;
bm(71).t1 = '2015-apr-01:07:29:24';
bm(71).t2 = '2015-apr-01:14:02:38';
bm(71).dk = 45;
bm(71).sched = '4_ffflat_mapo_40a_dk045_002.sch';
bm(71).filename = '20150401_072939';
bm(71).notes = '';

return


