
package modal_variations

import "core:math"
import "core:math/rand"
import "core:mem"
import "core:sort"
import "core:os"
import "core:fmt"

main :: proc() {
    using fmt;

    N :: FFT_SIZE * 2;
    data: [N]f32;

    get_sine :: proc(n: int, frac: f32) -> f32 {
        return math.sin(cast(f32) n * math.TAU / frac);
    }

    for n in 0..<N do data[n] = get_sine(n, 4.8) + get_sine(n, 8) + 0.6;

    modal_analysis(data[:], 24);

    println();
    println();
    println();

    noise: [cast(int) SAMPLE_RATE * 16]f32;
    for _, n in noise do noise[n] = rand.float32_range(-1, 1);
    println(noise[0:16]);

    apply_nonmodal_variation(noise[:], 5, Range{4, 12}, Range{2, 8});

    println(noise[0:16]);
    success := os.write_entire_file("Noise", transmute([]byte) noise[:]);
    println(success);
}


debug :: false;

SAMPLE_RATE : f32 : 44100;
FFT_SIZE    :     : 2048;



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

biquad_calculate_coefficients :: proc(frequency: f32, attenuation_db: f32, quality: f32) -> (coeffs: Biquad_Coefficients) {
    ω := math.TAU * frequency / SAMPLE_RATE;
    cos_ω := math.cos(ω);

    α := math.sin(ω) / (2 * quality);
    A := math.pow(10, -attenuation_db / 40);

    coeffs.b0 = 1  + α * A;
    coeffs.b1 = -2 * cos_ω;
    coeffs.b2 = 1  - α * A;

    coeffs.a0 = 1  + α / A;
    coeffs.a1 = coeffs.b1;
    coeffs.a2 = 1  - α / A;

    return;
}

Range :: struct {
    min, max: f32,
}

apply_nonmodal_variation :: proc(buffer: []f32, num_filters: int, gain: Range, q: Range) {
    for _ in 0..<num_filters {
        log_f := rand.float32_range(0, 1);

        frequency := (log_f * log_f) * 20_000;
        attenuation := rand.float32_range(gain.min, gain.max);
        quality     := rand.float32_range(q.min, q.max);

        coeffs := biquad_calculate_coefficients(frequency, attenuation, quality);

        state: Biquad_State;
        for value, n in buffer do buffer[n] = biquad_process(value, &state, coeffs);
    }
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
    for value, n in data do data[n] = complex(imag(value) / N, real(value) / N);
    cooley_tukey(data);
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

    when debug {
        println(magnitudes[:len(data)/2+1]);
        println();
    }

    for n, last_minimum_index: int; n < len(magnitudes) / 2; n += 1 {

        greater_than_previous: bool;
        if n == 0 do greater_than_previous = magnitudes[n] > 0;
        else      do greater_than_previous = magnitudes[n] > magnitudes[n-1];

        greater_than_next := magnitudes[n] > magnitudes[n+1];

        when debug do println(magnitudes[n], greater_than_previous, greater_than_next);

        if !greater_than_previous && !greater_than_next do last_minimum_index = n;
        else if greater_than_previous && greater_than_next {
            local_maximum := n;

            for greater_than_previous || greater_than_next {
                n += 1;
                greater_than_previous = magnitudes[n] > magnitudes[n-1];
                greater_than_next     = magnitudes[n] > magnitudes[n+1];
                when debug do println(magnitudes[n], greater_than_previous, greater_than_next);
            }
            // if n >= len(magnitudes) / 2 do break;

            last_minimum_value := magnitudes[last_minimum_index];
            if last_minimum_index == 0 do last_minimum_value = 0;

            when debug do println("new peak: ", last_minimum_value, " ", magnitudes[n]);

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
    envelope: [16]f32,
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

            copy(spectral_frames[l][:], alloc_for_r2c_fft(window[:]));
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
        if debug do println("\n\n\n");
        for all_peaks_in_one_frame in all_peaks_in_all_frames {
            println("\n=====================");
            for peak in all_peaks_in_one_frame do if peak.weight > 0.01 do println(peak);
        }
    }

    peaks_to_track := make([]Spectral_Peak, min(mode_count, len(all_peaks_in_all_frames[2])));
    copy(peaks_to_track, all_peaks_in_all_frames[2][:]);
    // for n in 0..<len(peaks) do peaks[n] = all_peaks[n];

    residue = make([]f32,  len(data));
    modes   = make([]Mode, len(peaks_to_track));

    for peak, n in peaks_to_track {

        bin := spectral_frames[0][peak.index];
        new_mode: Mode;
        new_mode.frequency = f32(peak.index) * SAMPLE_RATE / cast(f32) FFT_SIZE;
        new_mode.phase = math.atan2(imag(bin), real(bin));
        // new_mode.amplitude = peak.weight / cast(f32) FFT_SIZE;

        for frame, l in spectral_frames {
            bin = frame[peak.index];
            new_mode.envelope[l] = abs(bin) / cast(f32) FFT_SIZE;
            // println(peak.weight / cast(f32) FFT_SIZE);
        }

        modes[n] = new_mode;
    }
    //    return modes[:];
    println("\n\n\n");
    for mode in modes do println(mode, "\n");
    return;
}
