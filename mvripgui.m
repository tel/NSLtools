function x = mvripgui(action)
% MURIPGUI Callback function for mvrip_gui.m to generate moving ripples

% Auther: Taishih Chi (tschi@isr.umd.edu), NSL, UMD
% v1.00: 06-Jan-99


h = get(gcf,'children');
Tags = get(h, 'tag');
BW = str2num(get(h(strmatch('EditText01', Tags)), 'string'));% Bandwidth
Ra = get(h(strmatch('EditText02', Tags)), 'string');% Relative Amp.
Rt = get(h(strmatch('EditText03', Tags)), 'string');% Rate
Sl = get(h(strmatch('EditText04', Tags)), 'string');% Scale
Ph = get(h(strmatch('EditText05', Tags)), 'string');% Phase
Du = str2num(get(h(strmatch('EditText06', Tags)), 'string'));% Duration
SF = str2num(get(h(strmatch('EditText07', Tags)), 'string'));% Sampling Freq.
df = str2num(get(h(strmatch('EditText08', Tags)), 'string'));% Freq. Spacing
R0 = str2num(get(h(strmatch('EditText09', Tags)), 'string'));% Roll-off
Am = str2num(get(h(strmatch('EditText10', Tags)), 'string'));% Amp. Mod.
AU = 2-get(h(strmatch('PopupMenu2', Tags)), 'value');	% Amp. Unit	
CS = 2-get(h(strmatch('PopupMenu1', Tags)), 'value');	% Comp. Spacing

% get paras for the multimvripfft input

Paras = ['Ra';'Rt';'Sl';'Ph'];
[m n] = size(Paras);

% get the total number set of parameters
No = zeros(1, m);
for i = 1:m
	eval(['No(' num2str(i) ') = length(str2num(' Paras(i,:) '));']);
end
[tol_no, ind] = max(No);

for i = 1:m
	if (i==ind)
		temp = Paras(i,:);
		eval(['m_beg = findstr(' temp ',''['');']);
		eval(['m_end = findstr(' temp ','']'');']);
		eval(['comb_no = length(str2num(' temp '(m_beg:m_end)));']); 
		eval(['begin_id = length(str2num(' temp '(1:m_beg-1)))+1;']);
		eval([temp '=str2num(' temp ');']);	
	elseif (i~=ind & No(i)<tol_no)
		eval([Paras(i,:) '=str2num(' Paras(i,:) '(1))*ones(1, tol_no);']);
	elseif (i~=ind & No(i)==tol_no)
		eval([Paras(i,:) '=str2num(' Paras(i,:) ');']);		
	end
end

f0 = round(SF/(2^(BW+1))*1000);				% lowest freq.
paras = [8 8 -2 log2(SF/16)];
SF = round(SF/16*16384);
Ph = Ph/180*pi;
conds = [Du, f0, BW, SF, CS, df, R0, AU, Am];


% generate moving ripple

x = [];
total_no = tol_no-comb_no+(comb_no~=0);
for i = 1:total_no;
	if (i==begin_id)
		para = zeros(comb_no, 4);
		for k = 0:max(comb_no-1,0)
			para(k+1,:) = [Ra(i+k), Rt(i+k), Sl(i+k), Ph(i+k)];
		end
		i = i+k;
	else
		para = [Ra(i), Rt(i), Sl(i), Ph(i)];	
	end 
	xtmp = multimvripfft(para, conds);
	x = [x;xtmp;zeros(round(0.05*SF), 1)];
end

switch(action)
	case 'play',
		if (SF~=8192), x = resample(x,8192,SF); end;
		soundsc(x);
	case 'plot',
		subplot(h(strmatch('Axes1', Tags)));
		cla; axis([0 1 0 1]);
		text(0.5,0.5,'Please Wait...','fontw','bold',...
			'hori','center','fontsize', 24 );
		drawnow;
		y = wav2aud(x, paras);
		aud_plot(y, paras);
		set(gca, 'Tag', 'Axes1');
end
