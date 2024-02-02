function [data, auxData, metaData, txtData, weights] = mydata_Lampropholis_delicata 

% http://www.debtheory.org/wiki/index.php?title=Mydata_file#Metadata

%% set metaData
metaData.phylum     = 'Chordata'; 
metaData.class      = 'Reptilia'; 
metaData.order      = 'Squamata'; 
metaData.family     = 'Agamidae';
metaData.species    = 'Lampropholis_delicata'; 
metaData.species_en = 'dark-flecked garden sunskink'; 
metaData.T_typical  = C2K(32); % K, body temp
metaData.data_0     = {'ab'; 'tp'; 'am'; 'Lb'; 'Lp'; 'Li'; 'Wwb'; 'Wwp'; 'Wwi'; 'Ri'}; % zero-variate data labels: http://www.debtheory.org/wiki/index.php?title=Zero-variate_data
metaData.data_1     = {'t-Lw', 'L-Ww'}; % uni-variate data labels:  http://www.debtheory.org/wiki/index.php?title=Univariate_data 

metaData.COMPLETE = 2.4; % using criteria of LikaKear2011 http://www.debtheory.org/wiki/index.php?title=Completeness

metaData.author   = {'Kristoffer Wild'};    
metaData.date_subm = [2024 01 16];              
metaData.email    = {'kristofferw@unimelb.edu.au'};             
metaData.address  = {'School of Biological Sciences, The University of Melbourne, Parkville, VIC 3010'};   

%% set data
% zero-variate data

data.ab = 29.3;    units.ab = 'd';    label.ab = 'age at birth';             bibkey.ab = 'Current study';   
temp.ab = C2K(28);  units.temp.ab = 'K'; label.temp.ab = 'temperature'; comment.ab = 'Raw data from our study';

data.tp = 5*30;    units.tp = 'd';    label.tp = 'time since birth at puberty';           bibkey.tp = '';
temp.tp = C2K(29);  units.temp.tp = 'K'; label.temp.tp = 'temperature'; comment.tp = 'Noble DB and size from Forsman & Shine 1995'; 

data.am = 365*3;    units.am = 'd';    label.am = 'life span';                bibkey.am = 'Wilson & Swan, 2010';   
temp.am = C2K(25);  units.temp.am = 'K'; label.temp.am = 'temperature'; comment.am = '2-4 years'; 

data.Lb  = 1.78;   units.Lb  = 'cm';  label.Lb  = 'total length at birth';   bibkey.Lb  = 'Current study'; comment.Lb = '';  
data.Lp  = 3.2;   units.Lp  = 'cm';  label.Lp  = 'total length at puberty'; bibkey.Lp  = 'Noble DB and Forsman & Shine 1995'; comment.Lp = '';
data.Li  = 4.7;   units.Li  = 'cm';  label.Li  = 'ultimate total length';   bibkey.Li  = 'Kar et al., 2023'; comment.Li = '';

data.Wwb = 0.116;   units.Wwb = 'g';   label.Wwb = 'wet weight at birth';     bibkey.Wwb = 'Current study'; comment.Wwb = '';  
data.Wwp = 0.74;   units.Wwp = 'g';   label.Wwp = 'wet weight at puberty';   bibkey.Wwp = 'Noble DB and Forsman & Shine 1995';comment.Wwp = '';
data.Wwi = 2.9;   units.Wwi = 'g';   label.Wwi = 'ultimate wet weight';     bibkey.Wwi = 'Kar et al., 2023'; comment.Wwi = '';

data.Ri  = 0.04;   units.Ri  = '#/d'; label.Ri  = 'maximum reprod rate';     bibkey.Ri  = ''; comment.Ri = 'Noble DB';
temp.Ri = C2K(29);    units.temp.Ri = 'K'; label.temp.Ri = 'temperature';

% uni-variate data - 


%L-W data
LWdata = csvread('Raw_data/LW_data.csv',1,0);
data.LW=LWdata;
units.LW = {'cm', 'g'};  label.LW = {'Length', 'Wet weight'};  
bibkey.LW = 'Noble DB'; comment.LW = '';

%t-L data
tLdata = csvread('Raw_data/tL_data.csv',1,0);
data.tL = tLdata;
units.tL = {'d', 'cm'};  label.tL = {'time since birth', 'Length'};  
temp.tL = C2K(26);  % this is the constant temperature equivalent
units.temp.tL = 'K'; label.temp.tL = 'temperature';
bibkey.tL = 'Kar et al., 2023'; comment.tL = '';

%t-W data
tWdata = csvread('Raw_data/tW_data.csv', 1,0);
data.tW = tWdata;
units.tW = {'d', 'g'};  label.tW = {'time since birth', 'Wet weight'};  
temp.tW = C2K(26);  % this is the constant temperature equivalent
units.temp.tW = 'K'; label.temp.tW = 'temperature';
bibkey.tW = 'Kar et al., 2023_a'; comment.tW = '';

% O2 consumption 
% O2data= csvread('Raw_data/WJO_data.csv');
% data.WJO = O2data;
% units.WJO   = {'g', 'ml O2/min'};  label.WJO = {'wet weight','O2 consumption'};  
% temp.WJO    = C2K(28);  units.temp.WJO = 'K'; label.temp.WJO = 'temperature';
% bibkey.WJO = 'Kar et al., 2023_b';



%% set weights for all real data
weights = setweights(data, []);
weights.tp = .5 * weights.tp;
weights.Lp = .5 * weights.Lp;
weights.Wwp = .5 * weights.Wwp;
% weights.WJO = .5 * weights.WJO;

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);

%% pack auxData and txtData for output
auxData.temp = temp;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
%txtData.comment = comment;





%% References
bibkey = 'Wiki'; type = 'Misc'; bib = ...%
'howpublished = {\url{http://en.wikipedia.org/wiki/my_pet}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Kooy2010'; type = 'Book'; bib = [ ...  % used in setting of chemical parameters and pseudodata
'author = {Kooijman, S.A.L.M.}, ' ...
'year = {2010}, ' ...
'title  = {Dynamic Energy Budget theory for metabolic organisation}, ' ...
'publisher = {Cambridge Univ. Press, Cambridge}, ' ...
'pages = {Table 4.2 (page 150), 8.1 (page 300)}, ' ...
'howpublished = {\url{http://www.bio.vu.nl/thb/research/bib/Kooy2010.html}}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
% bibkey = 'Punzo2001'; type = 'Article'; bib = [ ... %
% 'author = {Punzo, F.}, ' ... 
% 'year = {2001}, ' ...
% 'title = {The Mediterranean Gecko, Hemidactylus turcicus: Life in an urban landscape}, ' ...
% 'journal = {Florida Scientist }, ' ...
% 'volume = {64}, ' ...
% 'pages = {56-66}'];
% metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
% %