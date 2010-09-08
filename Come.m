yposition = 10;

str = ['c  ';'o  ';'me ';'h  ';'o  ';'me ';'r  ';'i  ';'ght';'a  ';'w  ';'ay '];
pos = [50 100 150 280 360 440 530 620 710 800 850 910];

for k = 1:12,
	text('position', [pos(k), yposition], 'str', deblank(str(k, :)), ...
		'fonts', 16, 'fontn', 'times', 'fontw', 'b', 'ho', 'ce');
end;

