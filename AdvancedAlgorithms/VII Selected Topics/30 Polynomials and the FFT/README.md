/* 
 * Free FFT and convolution (Java)
 


public final class Fft {
	
	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This is a wrapper function.
	 */
	public static void transform(double[] real, double[] imag) {
		int n = real.length;
		if (n != imag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		if (n == 0)
			return;
		else if ((n & (n - 1)) == 0)  // Is power of 2
			transformRadix2(real, imag);
		else  // More complicated algorithm for arbitrary sizes
			transformBluestein(real, imag);
	}
	
	
	/* 
	 * Computes the inverse discrete Fourier transform (IDFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This is a wrapper function. This transform does not perform scaling, so the inverse is not a true inverse.
	 */
	public static void inverseTransform(double[] real, double[] imag) {
		transform(imag, real);
	}
	
	
	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector's length must be a power of 2. Uses the Cooley-Tukey decimation-in-time radix-2 algorithm.
	 */
	public static void transformRadix2(double[] real, double[] imag) {
		// Length variables
		int n = real.length;
		if (n != imag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		int levels = 31 - Integer.numberOfLeadingZeros(n);  // Equal to floor(log2(n))
		if (1 << levels != n)
			throw new IllegalArgumentException("Length is not a power of 2");
		
		// Trigonometric tables
		double[] cosTable = new double[n / 2];
		double[] sinTable = new double[n / 2];
		for (int i = 0; i < n / 2; i++) {
			cosTable[i] = Math.cos(2 * Math.PI * i / n);
			sinTable[i] = Math.sin(2 * Math.PI * i / n);
		}
		
		// Bit-reversed addressing permutation
		for (int i = 0; i < n; i++) {
			int j = Integer.reverse(i) >>> (32 - levels);
			if (j > i) {
				double temp = real[i];
				real[i] = real[j];
				real[j] = temp;
				temp = imag[i];
				imag[i] = imag[j];
				imag[j] = temp;
			}
		}
		
		// Cooley-Tukey decimation-in-time radix-2 FFT
		for (int size = 2; size <= n; size *= 2) {
			int halfsize = size / 2;
			int tablestep = n / size;
			for (int i = 0; i < n; i += size) {
				for (int j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
					int l = j + halfsize;
					double tpre =  real[l] * cosTable[k] + imag[l] * sinTable[k];
					double tpim = -real[l] * sinTable[k] + imag[l] * cosTable[k];
					real[l] = real[j] - tpre;
					imag[l] = imag[j] - tpim;
					real[j] += tpre;
					imag[j] += tpim;
				}
			}
			if (size == n)  // Prevent overflow in 'size *= 2'
				break;
		}
	}
	
	
	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This requires the convolution function, which in turn requires the radix-2 FFT function.
	 * Uses Bluestein's chirp z-transform algorithm.
	 */
	public static void transformBluestein(double[] real, double[] imag) {
		// Find a power-of-2 convolution length m such that m >= n * 2 + 1
		int n = real.length;
		if (n != imag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		if (n >= 0x20000000)
			throw new IllegalArgumentException("Array too large");
		int m = Integer.highestOneBit(n) * 4;
		
		// Trignometric tables
		double[] cosTable = new double[n];
		double[] sinTable = new double[n];
		for (int i = 0; i < n; i++) {
			int j = (int)((long)i * i % (n * 2));  // This is more accurate than j = i * i
			cosTable[i] = Math.cos(Math.PI * j / n);
			sinTable[i] = Math.sin(Math.PI * j / n);
		}
		
		// Temporary vectors and preprocessing
		double[] areal = new double[m];
		double[] aimag = new double[m];
		for (int i = 0; i < n; i++) {
			areal[i] =  real[i] * cosTable[i] + imag[i] * sinTable[i];
			aimag[i] = -real[i] * sinTable[i] + imag[i] * cosTable[i];
		}
		double[] breal = new double[m];
		double[] bimag = new double[m];
		breal[0] = cosTable[0];
		bimag[0] = sinTable[0];
		for (int i = 1; i < n; i++) {
			breal[i] = breal[m - i] = cosTable[i];
			bimag[i] = bimag[m - i] = sinTable[i];
		}
		
		// Convolution
		double[] creal = new double[m];
		double[] cimag = new double[m];
		convolve(areal, aimag, breal, bimag, creal, cimag);
		
		// Postprocessing
		for (int i = 0; i < n; i++) {
			real[i] =  creal[i] * cosTable[i] + cimag[i] * sinTable[i];
			imag[i] = -creal[i] * sinTable[i] + cimag[i] * cosTable[i];
		}
	}
	
	
	/* 
	 * Computes the circular convolution of the given real vectors. Each vector's length must be the same.
	 */
	public static void convolve(double[] x, double[] y, double[] out) {
		int n = x.length;
		if (n != y.length || n != out.length)
			throw new IllegalArgumentException("Mismatched lengths");
		convolve(x, new double[n], y, new double[n], out, new double[n]);
	}
	
	
	/* 
	 * Computes the circular convolution of the given complex vectors. Each vector's length must be the same.
	 */
	public static void convolve(double[] xreal, double[] ximag,
			double[] yreal, double[] yimag, double[] outreal, double[] outimag) {
		
		int n = xreal.length;
		if (n != ximag.length || n != yreal.length || n != yimag.length
				|| n != outreal.length || n != outimag.length)
			throw new IllegalArgumentException("Mismatched lengths");
		
		xreal = xreal.clone();
		ximag = ximag.clone();
		yreal = yreal.clone();
		yimag = yimag.clone();
		transform(xreal, ximag);
		transform(yreal, yimag);
		
		for (int i = 0; i < n; i++) {
			double temp = xreal[i] * yreal[i] - ximag[i] * yimag[i];
			ximag[i] = ximag[i] * yreal[i] + xreal[i] * yimag[i];
			xreal[i] = temp;
		}
		inverseTransform(xreal, ximag);
		
		for (int i = 0; i < n; i++) {  // Scaling (because this FFT implementation omits it)
			outreal[i] = xreal[i] / n;
			outimag[i] = ximag[i] / n;
		}
	}
	
}
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51
52
53
54
55
56
57
58
59
60
61
62
63
64
65
66
67
68
69
70
71
72
73
74
75
76
77
78
79
80
81
82
83
84
85
86
87
88
89
90
91
92
93
94
95
96
97
98
99
100
101
102
103
104
105
106
107
108
109
110
111
112
113
114
115
116
117
118
119
120
121
122
123
124
125
126
127
128
129
130
131
132
133
134
135
136
137
138
139
140
141
142
143
144
145
146
147
148
149
150
151
152
153
154
155
156
157
158
159
160
161
162
163
164
165
166
167
168
169
170
171
172
173
174
175
176
177
178
179
180
181
182
183
/* 
 * Free FFT and convolution (Java)
 
 
 
public final class Fft {
	
	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This is a wrapper function.
	 */

	public static void transform(double[] real, double[] imag) {

		int n = real.length;

		if (n != imag.length)

			throw new IllegalArgumentException("Mismatched lengths");

		if (n == 0)

			return;

		else if ((n & (n - 1)) == 0)  // Is power of 2

			transformRadix2(real, imag);

		else  // More complicated algorithm for arbitrary sizes

			transformBluestein(real, imag);

	}
	
	

	/* 
	 * Computes the inverse discrete Fourier transform (IDFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This is a wrapper function. This transform does not perform scaling, so the inverse is not a true inverse.
	 */

	public static void inverseTransform(double[] real, double[] imag) {

		transform(imag, real);

	}
	
	

	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector's length must be a power of 2. Uses the Cooley-Tukey decimation-in-time radix-2 algorithm.
	 */

	public static void transformRadix2(double[] real, double[] imag) {

		// Length variables

		int n = real.length;

		if (n != imag.length)

			throw new IllegalArgumentException("Mismatched lengths");

		int levels = 31 - Integer.numberOfLeadingZeros(n);  // Equal to floor(log2(n))

		if (1 << levels != n)

			throw new IllegalArgumentException("Length is not a power of 2");
		

		// Trigonometric tables

		double[] cosTable = new double[n / 2];

		double[] sinTable = new double[n / 2];

		for (int i = 0; i < n / 2; i++) {

			cosTable[i] = Math.cos(2 * Math.PI * i / n);

			sinTable[i] = Math.sin(2 * Math.PI * i / n);

		}
		

		// Bit-reversed addressing permutation

		for (int i = 0; i < n; i++) {

			int j = Integer.reverse(i) >>> (32 - levels);

			if (j > i) {

				double temp = real[i];

				real[i] = real[j];

				real[j] = temp;

				temp = imag[i];

				imag[i] = imag[j];

				imag[j] = temp;

			}

		}
		

		// Cooley-Tukey decimation-in-time radix-2 FFT

		for (int size = 2; size <= n; size *= 2) {

			int halfsize = size / 2;

			int tablestep = n / size;

			for (int i = 0; i < n; i += size) {

				for (int j = i, k = 0; j < i + halfsize; j++, k += tablestep) {

					int l = j + halfsize;

					double tpre =  real[l] * cosTable[k] + imag[l] * sinTable[k];

					double tpim = -real[l] * sinTable[k] + imag[l] * cosTable[k];

					real[l] = real[j] - tpre;

					imag[l] = imag[j] - tpim;

					real[j] += tpre;

					imag[j] += tpim;

				}

			}

			if (size == n)  // Prevent overflow in 'size *= 2'

				break;

		}

	}
	
	

	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This requires the convolution function, which in turn requires the radix-2 FFT function.
	 * Uses Bluestein's chirp z-transform algorithm.
	 */

	public static void transformBluestein(double[] real, double[] imag) {

		// Find a power-of-2 convolution length m such that m >= n * 2 + 1

		int n = real.length;

		if (n != imag.length)

			throw new IllegalArgumentException("Mismatched lengths");

		if (n >= 0x20000000)

			throw new IllegalArgumentException("Array too large");

		int m = Integer.highestOneBit(n) * 4;
		

		// Trignometric tables

		double[] cosTable = new double[n];

		double[] sinTable = new double[n];

		for (int i = 0; i < n; i++) {

			int j = (int)((long)i * i % (n * 2));  // This is more accurate than j = i * i

			cosTable[i] = Math.cos(Math.PI * j / n);

			sinTable[i] = Math.sin(Math.PI * j / n);

		}
		

		// Temporary vectors and preprocessing

		double[] areal = new double[m];

		double[] aimag = new double[m];

		for (int i = 0; i < n; i++) {

			areal[i] =  real[i] * cosTable[i] + imag[i] * sinTable[i];

			aimag[i] = -real[i] * sinTable[i] + imag[i] * cosTable[i];

		}

		double[] breal = new double[m];

		double[] bimag = new double[m];

		breal[0] = cosTable[0];

		bimag[0] = sinTable[0];

		for (int i = 1; i < n; i++) {

			breal[i] = breal[m - i] = cosTable[i];

			bimag[i] = bimag[m - i] = sinTable[i];

		}
		

		// Convolution

		double[] creal = new double[m];

		double[] cimag = new double[m];

		convolve(areal, aimag, breal, bimag, creal, cimag);
		

		// Postprocessing

		for (int i = 0; i < n; i++) {

			real[i] =  creal[i] * cosTable[i] + cimag[i] * sinTable[i];

			imag[i] = -creal[i] * sinTable[i] + cimag[i] * cosTable[i];

		}

	}
	
	

	/* 
	 * Computes the circular convolution of the given real vectors. Each vector's length must be the same.
	 */

	public static void convolve(double[] x, double[] y, double[] out) {

		int n = x.length;

		if (n != y.length || n != out.length)

			throw new IllegalArgumentException("Mismatched lengths");

		convolve(x, new double[n], y, new double[n], out, new double[n]);

	}
	
	

	/* 
	 * Computes the circular convolution of the given complex vectors. Each vector's length must be the same.
	 */

	public static void convolve(double[] xreal, double[] ximag,

			double[] yreal, double[] yimag, double[] outreal, double[] outimag) {
		

		int n = xreal.length;

		if (n != ximag.length || n != yreal.length || n != yimag.length

				|| n != outreal.length || n != outimag.length)

			throw new IllegalArgumentException("Mismatched lengths");
		

		xreal = xreal.clone();

		ximag = ximag.clone();

		yreal = yreal.clone();

		yimag = yimag.clone();

		transform(xreal, ximag);

		transform(yreal, yimag);
		

		for (int i = 0; i < n; i++) {

			double temp = xreal[i] * yreal[i] - ximag[i] * yimag[i];

			ximag[i] = ximag[i] * yreal[i] + xreal[i] * yimag[i];

			xreal[i] = temp;

		}

		inverseTransform(xreal, ximag);
		

		for (int i = 0; i < n; i++) {  // Scaling (because this FFT implementation omits it)

			outreal[i] = xreal[i] / n;

			outimag[i] = ximag[i] / n;

		}

	}
	
}
تبدیل فوریه سریع در زبان متلب
با استفاده از برنامه زیر می‌توان تبدیل فوریه سریع را در زبان متلب محاسبه کرد.

MATLAB
function [spectrum, freq] = autofft(xs, ts, fftset)

%% nargin check
if nargin < 2
    error('Not enough input arguments.');
elseif nargin > 3
    error('Too many input arguments.');
end
%
%% Convert row vectors to column vectors if needed
if size(xs, 1) == 1         % samples
    xs = xs(:);                     
end
if size(ts(:), 1) == 1      % sampling frequency
	fs = ts;                    
else
    fs = 1 / (ts(2) - ts(1)); 
end
%
%% Specify default fftset
defset = struct('nwin', size(xs, 1), ...
                'twin', size(xs, 1) / fs, ...
                'overlap', 50, ...
                'lowpass', fs/2, ...
                'window', 'u', ...
                'averaging', 'lin', ...
                'jw', '1', ...
                'unit', 'pow');
deffields = fieldnames(defset);
%
%% Set analyser parameters
if nargin == 2  % use default fftset
    fftset = defset;
else            % use user-defined fftset  
    % Check whether there is user-defined 'nwin' or 'twin' parameter
    if isfield(fftset, 'nwin')
        fftset.twin = fftset.nwin / fs;
    elseif isfield(fftset, 'twin')
        fftset.nwin = round(fftset.twin * fs);
    end
    % Set unspecified parameters to default
    for i = 1:numel(deffields)        
        if isfield(fftset, deffields{i}) == 0
            fftset.(deffields{i}) = defset.(deffields{i});
        end
    end
end
% Generate frequency vector
freq = (fs * (0:(fftset.nwin/2)) / fftset.nwin)';
% Set allowed frequencies for the low-pass filtering
freq = freq(freq <= fftset.lowpass);
% Set handling of the last spectral line
if freq(end) == fs/2
    sh = 1; % Nyquist frequency is at the last spectral line
else
	sh = 0; % Nyquist frequency is not at the last spectral line
end
% Set number of overlaping samples
fftset.overlap = round(fftset.nwin * fftset.overlap / 100);
% Set indices for the signal segmentation
imax = floor((size(xs, 1)-fftset.overlap) / (fftset.nwin-fftset.overlap));
                                                % number of windows
ind = zeros(imax,2);                            % matrix of indices
ni = 1;                                         % pointer
for i = 1:imax                                  % cycle through windows
    ind(i,1) = ni; 
    ni = ni + fftset.nwin - 1;
    ind(i,2) = ni;
    ni = ni - fftset.overlap + 1;
end
% Generate the time weighting window
[fftset.window, fftset.noiseband] = timeweight(fftset.window, fftset.nwin, fs);
% Set constant for the jw weigthing
switch lower(fftset.jw)
    case '1/jw2'
        fftset.jw = - 1 ./ (4 * pi^2 * freq.^2);
    case '1/jw'
        fftset.jw = 1 ./ (2i * pi * freq);
    case 'jw'
        fftset.jw = 2i * pi * freq;
    case 'jw2'
        fftset.jw = - 4 * pi^2 * freq.^2;
    otherwise
        fftset.jw = 1;
end
%
%% Spectral analyser
spectrum = zeros(size(freq, 1), size(xs, 2));
%
for i = 1:size(xs, 2)
    % Preallocate an array for temporary spectra
    spectra = zeros(size(freq, 1), imax);
    % Frequency analysis of individual segments
    for j = 1:imax
        % Fast Fourier transformation of the weighted segment
        fftSamp = fft(fftset.window .* xs(ind(j,1):ind(j,2), i), ...
                      fftset.nwin) / fftset.nwin;
        % Application of the jw weigthing and the low-pass filtering
        spectra(:, j) = fftset.jw .* fftSamp(1:size(freq, 1));
        % Evaluation of spectral unit
        switch lower(fftset.unit)
            case 'rms'           % Linear spectrum with rms magnitude
                spectra(2:end-sh,j) = (2/sqrt(2))*abs(spectra(2:end-sh,j));
                spectra(end-sh+1:end,j) = abs(spectra(end-sh+1:end,j))/sqrt(2);
            case 'pk'            % Linear spectrum with 0-peak magnitude
                spectra(2:end-sh,j) = 2 * abs(spectra(2:end-sh,j));
            case 'pp'            % Linear spectrum with peak-peak magnitude
                spectra(2:end-sh,j) = 4 * abs(spectra(2:end-sh,j));
                spectra(end-sh+1:end,j) = 2 * abs(spectra(end-sh+1:end,j));
            case 'psd'           % Power spectral density
                spectra(2:end-sh,j) = sqrt(2) * spectra(2:end-sh,j);
                spectra(:,j)        = (1 / fftset.noiseband) * spectra(:,j) ...
                                       .* conj(spectra(:,j));
            case {'rsd','rmspsd'}% Root mean square of PSD 
                spectra(2:end-sh,j) = sqrt(2) * spectra(2:end-sh,j);
                spectra(:,j)        = sqrt((1 / fftset.noiseband) * ...
                                       spectra(:,j) .* conj(spectra(:,j)));                
            otherwise            % Autospectrum
                spectra(2:end-sh,j) = sqrt(2) * spectra(2:end-sh,j);
                spectra(:,j)        = spectra(:,j) .* conj(spectra(:,j));
        end
    end
    % Spectral averaging
    switch lower(fftset.averaging)
        case {'energy', 'rms'}   % Energy averaging
            spectrum(:, i) = sqrt(sum(spectra.^2, 2) ./ imax);
        case {'max', 'peak'}     % Maximum peak hold averaging
            spectrum(:, i) = max(spectra, [], 2);
        case 'min'              % Minimum peak hold averaging
            spectrum(:, i) = min(spectra, [], 2);
        case 'none'             % No averaging
            if size(xs, 2) == 1
                spectrum = spectra;
            else
                spectrum{i} = spectra;
            end
        otherwise               % Linear averaging
            spectrum(:, i) = mean(spectra, 2);
    end
    % Remove imaginary part (residue due to spectra .* conj(spectra))
    spectrum = real(spectrum);
end
% End of main fucntion
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunction weightsig generates the time weighting window
%
% Input
%   - sym - symbol for the time weighting window
%   - n   - length of the time weighting window (samples)
%   - fs  - sampling frequency (Hz) 
%
% Output
%   - window    - the time weighting window
%   - noiseband - the noise power bandwidth of the time weighting window
%  
function [window, noiseband] = timeweight(sym, n, fs)
    % Generate the specified time weighting window
    switch sym(1)
        case 'b'    % Blackmann-Harris
            window = blackmanharris(n);
        case 'f'    % flat-top
            window = flattopwin(n);
        case 'h'    % Hann
            window = hann(n);
        case 'k'    % Kaiser-Bessel
            if length(sym) == 1
                % window with default beta = 0.5
                window = kaiser(n, 0.5);
            else
                % window with user specified beta
                window = kaiser(n, str2double(sym(2:end)));
            end
        case 'm'    % Hamming
            window = hamming(n);
        otherwise   % uniform
            window = rectwin(n);
    end
    % Adjust magnitude of the generated time weighting window
    window = window / mean(window);
    % Calculate the noise power bandwidth in Hz
    noiseband = enbw(window, fs);
end
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51
52
53
54
55
56
57
58
59
60
61
62
63
64
65
66
67
68
69
70
71
72
73
74
75
76
77
78
79
80
81
82
83
84
85
86
87
88
89
90
91
92
93
94
95
96
97
98
99
100
101
102
103
104
105
106
107
108
109
110
111
112
113
114
115
116
117
118
119
120
121
122
123
124
125
126
127
128
129
130
131
132
133
134
135
136
137
138
139
140
141
142
143
144
145
146
147
148
149
150
151
152
153
154
155
156
157
158
159
160
161
162
163
164
165
166
167
168
169
170
171
172
173
174
175
176
177
178
179
180
181
182
183
184
185
186
187
188
189
190
191
192
193
194
195
196
197
198
199
200
201
202
203
204
205
206
207
208
209
210
211
212
213
214
215

function [spectrum, freq] = autofft(xs, ts, fftset)
 

%% nargin check
 

if nargin < 2

    error('Not enough input arguments.');

elseif nargin > 3

    error('Too many input arguments.');
end
%
 

%% Convert row vectors to column vectors if needed
 

if size(xs, 1) == 1         % samples

    xs = xs(:);                     
end

if size(ts(:), 1) == 1      % sampling frequency

	fs = ts;                    
else

    fs = 1 / (ts(2) - ts(1)); 
end
%
 

%% Specify default fftset
 

defset = struct('nwin', size(xs, 1), ...

                'twin', size(xs, 1) / fs, ...

                'overlap', 50, ...

                'lowpass', fs/2, ...

                'window', 'u', ...

                'averaging', 'lin', ...

                'jw', '1', ...

                'unit', 'pow');

deffields = fieldnames(defset);
%
 

%% Set analyser parameters
 

if nargin == 2  % use default fftset

    fftset = defset;

else            % use user-defined fftset  

    % Check whether there is user-defined 'nwin' or 'twin' parameter

    if isfield(fftset, 'nwin')

        fftset.twin = fftset.nwin / fs;

    elseif isfield(fftset, 'twin')

        fftset.nwin = round(fftset.twin * fs);

    end

    % Set unspecified parameters to default

    for i = 1:numel(deffields)        

        if isfield(fftset, deffields{i}) == 0

            fftset.(deffields{i}) = defset.(deffields{i});

        end

    end
end

% Generate frequency vector
 

freq = (fs * (0:(fftset.nwin/2)) / fftset.nwin)';
% Set allowed frequencies for the low-pass filtering
freq = freq(freq <= fftset.lowpass);
% Set handling of the last spectral line
if freq(end) == fs/2
    sh = 1; % Nyquist frequency is at the last spectral line
else
	sh = 0; % Nyquist frequency is not at the last spectral line
end
% Set number of overlaping samples
fftset.overlap = round(fftset.nwin * fftset.overlap / 100);
% Set indices for the signal segmentation
imax = floor((size(xs, 1)-fftset.overlap) / (fftset.nwin-fftset.overlap));
                                                % number of windows
ind = zeros(imax,2);                            % matrix of indices
ni = 1;                                         % pointer
for i = 1:imax                                  % cycle through windows
    ind(i,1) = ni; 
    ni = ni + fftset.nwin - 1;
    ind(i,2) = ni;
    ni = ni - fftset.overlap + 1;
end
% Generate the time weighting window
[fftset.window, fftset.noiseband] = timeweight(fftset.window, fftset.nwin, fs);
% Set constant for the jw weigthing
switch lower(fftset.jw)

    case '1/jw2'
        fftset.jw = - 1 ./ (4 * pi^2 * freq.^2);

    case '1/jw'
        fftset.jw = 1 ./ (2i * pi * freq);

    case 'jw'
        fftset.jw = 2i * pi * freq;

    case 'jw2'
        fftset.jw = - 4 * pi^2 * freq.^2;
    otherwise
        fftset.jw = 1;
end
%
%% Spectral analyser
spectrum = zeros(size(freq, 1), size(xs, 2));
%
for i = 1:size(xs, 2)
    % Preallocate an array for temporary spectra
    spectra = zeros(size(freq, 1), imax);
    % Frequency analysis of individual segments
    for j = 1:imax
        % Fast Fourier transformation of the weighted segment
        fftSamp = fft(fftset.window .* xs(ind(j,1):ind(j,2), i), ...
                      fftset.nwin) / fftset.nwin;
        % Application of the jw weigthing and the low-pass filtering
        spectra(:, j) = fftset.jw .* fftSamp(1:size(freq, 1));
        % Evaluation of spectral unit
        switch lower(fftset.unit)

            case 'rms'           % Linear spectrum with rms magnitude
                spectra(2:end-sh,j) = (2/sqrt(2))*abs(spectra(2:end-sh,j));
                spectra(end-sh+1:end,j) = abs(spectra(end-sh+1:end,j))/sqrt(2);

            case 'pk'            % Linear spectrum with 0-peak magnitude
                spectra(2:end-sh,j) = 2 * abs(spectra(2:end-sh,j));

            case 'pp'            % Linear spectrum with peak-peak magnitude
                spectra(2:end-sh,j) = 4 * abs(spectra(2:end-sh,j));
                spectra(end-sh+1:end,j) = 2 * abs(spectra(end-sh+1:end,j));

            case 'psd'           % Power spectral density
                spectra(2:end-sh,j) = sqrt(2) * spectra(2:end-sh,j);
                spectra(:,j)        = (1 / fftset.noiseband) * spectra(:,j) ...
                                       .* conj(spectra(:,j));

            case {'rsd','rmspsd'}% Root mean square of PSD 
                spectra(2:end-sh,j) = sqrt(2) * spectra(2:end-sh,j);
                spectra(:,j)        = sqrt((1 / fftset.noiseband) * ...
                                       spectra(:,j) .* conj(spectra(:,j)));                
            otherwise            % Autospectrum
                spectra(2:end-sh,j) = sqrt(2) * spectra(2:end-sh,j);
                spectra(:,j)        = spectra(:,j) .* conj(spectra(:,j));
        end
    end
    % Spectral averaging
    switch lower(fftset.averaging)

        case {'energy', 'rms'}   % Energy averaging
            spectrum(:, i) = sqrt(sum(spectra.^2, 2) ./ imax);

        case {'max', 'peak'}     % Maximum peak hold averaging
            spectrum(:, i) = max(spectra, [], 2);

        case 'min'              % Minimum peak hold averaging
            spectrum(:, i) = min(spectra, [], 2);

        case 'none'             % No averaging
            if size(xs, 2) == 1
                spectrum = spectra;
            else
                spectrum{i} = spectra;
            end
        otherwise               % Linear averaging
            spectrum(:, i) = mean(spectra, 2);
    end
    % Remove imaginary part (residue due to spectra .* conj(spectra))
    spectrum = real(spectrum);
end
% End of main fucntion
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunction weightsig generates the time weighting window
%
% Input
%   - sym - symbol for the time weighting window
%   - n   - length of the time weighting window (samples)
%   - fs  - sampling frequency (Hz) 
%
% Output
%   - window    - the time weighting window
%   - noiseband - the noise power bandwidth of the time weighting window
%  
function [window, noiseband] = timeweight(sym, n, fs)
    % Generate the specified time weighting window
    switch sym(1)

        case 'b'    % Blackmann-Harris
            window = blackmanharris(n);

        case 'f'    % flat-top
            window = flattopwin(n);

        case 'h'    % Hann
            window = hann(n);

        case 'k'    % Kaiser-Bessel
            if length(sym) == 1
                % window with default beta = 0.5
                window = kaiser(n, 0.5);
            else
                % window with user specified beta
                window = kaiser(n, str2double(sym(2:end)));
            end

        case 'm'    % Hamming

            window = hamming(n);

        otherwise   % uniform

            window = rectwin(n);

    end

    % Adjust magnitude of the generated time weighting window

    window = window / mean(window);

    % Calculate the noise power bandwidth in Hz

    noiseband = enbw(window, fs);
end

