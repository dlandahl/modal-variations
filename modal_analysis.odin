
package modal_variations

import "core:math"
import "core:math/rand"
import "core:math/bits"
import "core:mem"
import "core:sort"
import "core:os"
import "core:fmt"

divisor :: inline proc() {
    fmt.println("\n\n==========================");
}

main :: proc() {
    using fmt;

    // noise: [cast(int) SAMPLE_RATE * 64]f32;
    // for _, n in noise do noise[n] = rand.float32_range(-1, 1);
    // 
    // apply_nonmodal_variation(noise[:], 12, Range{2, 8}, Range{2,6});
    // 
    // buffer := wave_encode(noise[:]);
    // defer delete(buffer);
    // 
    // os.write_entire_file("Test file.wav", transmute([]u8) buffer);

    // buffer, error := os.read_entire_file("generativei24.wav");
    // raw_data, header, success := wave_decode(buffer);

    raw_data: [1 << 18]f32;
    for value, n in raw_data {
        raw_data[n] = math.sin(cast(f32) n * math.TAU * 440 / cast(f32) SAMPLE_RATE);
    }

    complex_data := alloc_for_r2c_fft(raw_data[:]);
    cooley_tukey(complex_data[:]);
    inverse_cooley_tukey(complex_data[:]);

    for value, n in raw_data do raw_data[n] = real(complex_data[n]);

    new_buffer := wave_encode(raw_data[:], 1);
    os.write_entire_file("FILE", transmute([]u8) raw_data[:]);

    // delete(residue);
    // delete(modes);
    // delete(new_buffer);
    // delete(raw_data);
    // delete(buffer);
}

DEBUG       :: false;
SAMPLE_RATE :: 44100;
FFT_SIZE    :: 2048;



Biquad_State :: struct {
    x1, x2,
    y1, y2: f32,
}

Biquad_Coefficients :: struct {
    a0, a1, a2,
    b0, b1, b2: f32,
}

biquad_process :: proc(x0: f32, using state: ^Biquad_State, using coeffs: Biquad_Coefficients) -> (y0: f32) {
    y0 = (b0*x0 + b1*x1 + b2*x2 - a1*y1 - a2*y2) / a0;

    x2 = x1;
    x1 = x0;

    y2 = y1;
    y1 = y0;
    return;
}

biquad_calculate_coefficients :: proc(frequency, attenuation_db, quality: f32) -> (coeffs: Biquad_Coefficients) {
    using math;
    ω := τ * frequency / SAMPLE_RATE;

    α := sin(ω) / (2 * quality);
    A := pow(10, -attenuation_db / 40);

    coeffs.b0 = 1  + α * A;
    coeffs.b1 = -2 * cos(ω);
    coeffs.b2 = 1  - α * A;

    coeffs.a0 = 1  + α / A;
    coeffs.a1 = coeffs.b1;
    coeffs.a2 = 1  - α / A;

    return;
}

Range :: struct {
    min, max: f32,
}

apply_nonmodal_variation :: proc(buffer: []f32, num_filters: int, gain, q: Range) {
    for _ in 0..<num_filters {
        log_f := rand.float32_range(0, 1);

        frequency := math.pow(1, log_f) * 20_000;
        attenuation := rand.float32_range(gain.min, gain.max);
        quality     := rand.float32_range(q.min, q.max);

        coeffs := biquad_calculate_coefficients(frequency, attenuation, quality);

        state: Biquad_State;
        for value, n in buffer do buffer[n] = biquad_process(value, &state, coeffs);
    }
}



WAVE_HEADER_SIZE :: 44;
#assert (WAVE_HEADER_SIZE == size_of(Wave_Header));

Wave_Format :: enum u16 {
    INTEGER = 1,
    FLOAT   = 3, 
}

Wave_Header :: struct #packed {
    chunk_id        : [4]byte,
    chunk_size      : u32,
    format          : [4]byte,
    subchunk_1_id   : [4]byte,
    subchunk_1_size : u32,
    audio_format    : Wave_Format,
    channels        : u16,
    sample_rate     : u32,
    byte_rate       : u32,
    block_align     : u16,
    bits_per_sample : u16,
    subchunk_2_id   : [4]byte,
    subchunk_2_size : u32,
}

wave_generate_header :: proc(num_samples: int, num_channels: int = 1) -> Wave_Header {
    using header: Wave_Header;
    bytes_per_sample := size_of(f32);

    chunk_id        = {'R', 'I', 'F', 'F'};
    chunk_size      = u32(36 + num_samples * bytes_per_sample);
    format          = {'W', 'A', 'V', 'E'};
    subchunk_1_id   = {'f', 'm', 't', ' '};
    subchunk_1_size = 16;
    audio_format    = Wave_Format.FLOAT;
    channels        = u16(num_channels);
    sample_rate     = u32(SAMPLE_RATE);
    byte_rate       = u32(SAMPLE_RATE * bytes_per_sample * num_channels);
    block_align     = u16(bytes_per_sample * num_channels);
    bits_per_sample = u16(bytes_per_sample * 8);
    subchunk_2_id   = {'d', 'a', 't', 'a'};
    subchunk_2_size = u32(num_samples * bytes_per_sample);

    return header;
}

wave_read_header :: proc(file: []byte) -> (header: Wave_Header) {
    header_bytes: [WAVE_HEADER_SIZE]byte;
    copy(header_bytes[:], file);
    header = transmute(Wave_Header) header_bytes;
    return;
}

wave_encode :: proc(buffer: []f32, num_channels: int = 1) -> (file: []byte) {
    file = make([]byte, WAVE_HEADER_SIZE + len(buffer) * size_of(f32));

    header := transmute([WAVE_HEADER_SIZE]byte) wave_generate_header(len(buffer), num_channels);
    copy(file, header[:]);
    copy(file[WAVE_HEADER_SIZE:], transmute([]byte) buffer);

    return;
}

wave_decode :: proc(file: []byte) -> (data: []f32, header: Wave_Header, success: bool) {

    header = wave_read_header(file);
    success = true;

    if     header.chunk_id      != {'R', 'I', 'F', 'F'}
        || header.format        != {'W', 'A', 'V', 'E'}
        || header.subchunk_1_id != {'f', 'm', 't', ' '}
        || header.subchunk_2_id != {'d', 'a', 't', 'a'} {
            fmt.println("Invalid wave header");
            success = false;
            return;
    }

    data_chunk := file[WAVE_HEADER_SIZE:];
    data = make([]f32, len(data_chunk) / int(header.bits_per_sample / 8));

    switch header.audio_format {
    case .FLOAT:
        copy(data, transmute([]f32) data_chunk);

    case .INTEGER:
        max_value := 1 << (header.bits_per_sample - 1) - 1;

        if header.bits_per_sample == 16 do for _, n in data {
            data[n] = cast(f32) (transmute([]i16) data_chunk)[n];
        }
        else if header.bits_per_sample == 24 {
            array := data_chunk[:];
            for n := 0; len(array) >= 3; n += 1 {
                data[n] = cast(f32) (u32(array[2]) << 16 | u32(array[1]) << 8 | u32(array[0]));
                array = array[3:];
            }
        }

        for _, n in data {
            data[n] /= f32(max_value);
            if data[n] >= 1 do data[n] -= 2;
        }
    case:
        fmt.println("Invalid wave header");
        success = false;
        return;
    }
    return;
}



cooley_tukey :: proc(data: []complex64) {

    N := len(data);
    assert(math.is_power_of_two(N), "FFT size must be a power of two");
    if N < 2 do return;

    M := N / 2;

    temp := make([]complex64, M);
    for n in 0..<M do temp[n]     = data[n * 2 + 1];
    for n in 0..<M do data[n]     = data[n * 2];
    for n in 0..<M do data[n + M] = temp[n];
    delete(temp);

    cooley_tukey(data[:M]);
    cooley_tukey(data[M:]);

    for k in 0..<M {
        even := data[k];
        odd  := data[k + M];

        theta := -math.TAU * f32(k) / f32(N);
        w     := complex(math.cos(theta), math.sin(theta)) * odd;

        data[k]     = even + w;
        data[k + M] = even - w;
    }
}

inverse_cooley_tukey :: proc(data: []complex64) {

    N := cast(f32) len(data);
    for value, n in data do data[n] = complex(imag(value), real(value));
    cooley_tukey(data);
    for value, n in data do data[n] = complex(imag(value) / N, real(value) / N);
}

alloc_for_r2c_fft :: proc(data: []f32) -> []complex64 {

    cdata := make([]complex64, len(data));
    for _, n in data do cdata[n] = complex(data[n], 0);
    return cdata;
}

apply_hamming_window :: proc(data: []$T) {
    a_0 :: 0.53836;
    a_1 :: 1 - a_0;

    N := cast(f32) len(data);
    for _, n in data {
        data[n] *= a_0 - a_1 * math.cos(math.TAU * f32(n) / N);
    }
}



Spectral_Peak :: struct {
    weight: f32,
    index: int,
}

find_and_alloc_spectral_peaks :: proc(data: []complex64) -> (peaks: [dynamic]Spectral_Peak) {
    using fmt;

    magnitudes := make([]f32, len(data));
    defer delete(magnitudes);

    for value, n in data do magnitudes[n] = abs(value);

    when DEBUG {
        println(magnitudes[:len(data)/2+1]);
        println();
    }

    for n, last_minimum_index: int; n < len(magnitudes) / 2; n += 1 {

        greater_than_previous: bool;
        if n == 0 do greater_than_previous = magnitudes[n] > 0;
        else      do greater_than_previous = magnitudes[n] > magnitudes[n-1];

        greater_than_next := magnitudes[n] > magnitudes[n+1];

        when DEBUG do println(magnitudes[n], greater_than_previous, greater_than_next);

        if !greater_than_previous && !greater_than_next do last_minimum_index = n;
        else if greater_than_previous && greater_than_next {
            local_maximum := n;

            for greater_than_previous || greater_than_next {
                n += 1;
                greater_than_previous = magnitudes[n] > magnitudes[n-1];
                greater_than_next     = magnitudes[n] > magnitudes[n+1];
                when DEBUG do println(magnitudes[n], greater_than_previous, greater_than_next);
            }

            last_minimum_value := magnitudes[last_minimum_index];
            if last_minimum_index == 0 do last_minimum_value = 0;

            when DEBUG do println("new peak: ", last_minimum_value, " ", magnitudes[n]);

            average_minimum := (last_minimum_value + magnitudes[n]) / 2;

            new_peak: Spectral_Peak;
            new_peak.index = local_maximum;
            new_peak.weight = magnitudes[local_maximum] - average_minimum;
            if (new_peak.weight / cast(f32) len(data) > 0.01) do append(&peaks, new_peak);

            last_minimum_index = n;
        }
    }

    return;
}

Mode :: struct {
    envelope: [1024]f32,
    frequency: f32,
    phase: f32,
}

modal_analysis :: proc(data: []f32, mode_count: int) -> (modes: []Mode, residue: []f32) {
    using fmt;


    M :: FFT_SIZE;
    HOP_SIZE :: FFT_SIZE / 4;

    num_spectral_frames := (len(data) + 1) / HOP_SIZE;
    spectral_frames := make([][M]complex64, num_spectral_frames);

    {
        padded_data := make([]f32, len(data) + M);
        defer delete(padded_data);
        copy(padded_data, data);

        for l in 0..<num_spectral_frames {
            window: [M]f32;
            index := l * HOP_SIZE;

            copy(window[:], padded_data[index:index + M]);
            apply_hamming_window(window[:]);

            complex_buffer := alloc_for_r2c_fft(window[:]);
            copy(spectral_frames[l][:], complex_buffer);
            delete(complex_buffer);
            cooley_tukey(spectral_frames[l][:]);
        }
    }

    all_peaks_in_all_frames := make([][dynamic]Spectral_Peak, num_spectral_frames);
    defer {
        for peaks in all_peaks_in_all_frames do delete(peaks);
        delete(all_peaks_in_all_frames);
    }

    for frame, n in spectral_frames {
        the_frame := frame;
        all_peaks_in_all_frames[n] = find_and_alloc_spectral_peaks(the_frame[:]);
    }

    {
        comparison_function :: proc(a, b: Spectral_Peak) -> int {
            return sort.compare_f32s(b.weight, a.weight);
        }

        for all_peaks_in_one_frame in all_peaks_in_all_frames {
            sort.bubble_sort_proc(all_peaks_in_one_frame[:], comparison_function);
        }
    }

    {
        when DEBUG do println("\n\n\n");
        for all_peaks_in_one_frame in all_peaks_in_all_frames {
            println("\n=====================");
            for peak in all_peaks_in_one_frame do if peak.weight > 0.01 do println(peak);
        }
    }

    peaks_to_track := make([]Spectral_Peak, min(mode_count, len(all_peaks_in_all_frames[2])));
    copy(peaks_to_track, all_peaks_in_all_frames[2][:]);

    residue = make([]f32,  len(data));
    modes   = make([]Mode, len(peaks_to_track));

    for peak, n in peaks_to_track {

        bin := spectral_frames[2][peak.index];
        new_mode: Mode;
        new_mode.frequency = f32(peak.index) * SAMPLE_RATE / FFT_SIZE;
        new_mode.phase = math.atan2(imag(bin), real(bin));

        for frame, l in spectral_frames {
            bin = frame[peak.index];
            new_mode.envelope[l] = abs(bin) / FFT_SIZE;
        }

        modes[n] = new_mode;
    }

    residue = make([]f32, len(data) + M);

    for frame, l in spectral_frames {
        index := l * HOP_SIZE;
        slice := residue[index:index + M];

        the_frame := frame;
        inverse_cooley_tukey(the_frame[:]);
        for n in 0..<M {
            slice[n] += real(frame[n]);
        }
    }

    when DEBUG {
        println("\n\n\n");
        for mode in modes do println(mode, "\n");
    }

    return;
}
