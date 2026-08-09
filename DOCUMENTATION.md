# lucky-star — dokumentacja

`lucky-star` analizuje zdjęcia astronomiczne zapisane w plikach FITS. Oblicza metryki jakości, może wygenerować wykresy i mapę jakości, a także filtrować oraz porządkować pliki według wybranej metryki.

## 1. Instalacja i uruchomienie

Projekt jest aplikacją Rust. Do kompilacji potrzebny jest zainstalowany Rust wraz z Cargo.

```bash
cargo build --release
```

Po kompilacji program znajduje się w `target/release/lucky-star` (na Windows: `target/release/lucky-star.exe`). W trakcie pracy nad projektem można uruchamiać wersję developerską:

```bash
cargo run -- [OPCJE] [ŚCIEŻKA]
```

Przykład uruchomienia gotowego programu:

```bash
./target/release/lucky-star ./sesja
```

Domyślna ścieżka to bieżący katalog (`.`). Program przetwarza pliki z rozszerzeniem `.fits` i `.FITS` znajdujące się bezpośrednio w podanym katalogu.

Pełną listę opcji można wyświetlić poleceniem:

```bash
lucky-star --help
```

## 2. Podstawowa analiza

Analiza katalogu:

```bash
lucky-star ./sesja
```

Analiza pojedynczego pliku FITS:

```bash
lucky-star ./sesja/klatka_0001.fits
```

Dla katalogu program analizuje obrazy, wypisuje statystyki metryk i — przy pełnym przetwarzaniu katalogu — zapisuje mapę jakości oraz wykresy metryk. Wyniki graficzne są zapisywane w analizowanym katalogu.

Opcja `--verbose` włącza dodatkowe informacje diagnostyczne:

```bash
lucky-star ./sesja --verbose
```

## 3. Metryki

Program obsługuje następujące metryki:

| Metryka | Znaczenie | Lepsza wartość |
|---|---|---|
| `quality` | ogólna jakość obrazu | większa |
| `fwhm` | szerokość profilu gwiazdy; mniejsza oznacza ostrzejszy obraz | mniejsza |
| `quality_star_pattern` | jakość dopasowania zdefiniowanego wzorca gwiazd | większa |
| `background` | poziom tła po korekcji | mniejsza |
| `star_brightness` | jasność zarejestrowanych gwiazd wzorca | większa |
| `snr` | stosunek sygnału gwiazd wzorca do szumu tła | większa |

Wartości zależne od wzorca (`quality_star_pattern`, `star_brightness`, `snr`) wymagają poprawnie wczytanego wzorca gwiazd.

## 4. Filtrowanie plików — `--filter`

Opcja `--filter` uruchamia filtrowanie i przenosi odrzucone pliki do nowego folderu `removed_*`. Pliki, które spełniają wszystkie warunki, pozostają w katalogu źródłowym.

### 4.1. Progi względne

Bez sufiksu `-absolute` próg jest mnożnikiem mediany danej metryki. Na przykład:

```bash
lucky-star ./sesja --filter --snr 0.8
```

zachowuje obrazy z SNR co najmniej równym `0.8 × mediana(SNR)`.

Dla metryk, dla których mniejsza wartość jest lepsza, filtr zachowuje wartości nie większe od obliczonego progu. Przykład:

```bash
lucky-star ./sesja --filter --fwhm 1.2
```

zachowuje obrazy o FWHM nie większym niż `1.2 × mediana(FWHM)`.

Można zastosować wiele filtrów jednocześnie. Obraz musi spełnić wszystkie warunki:

```bash
lucky-star ./sesja --filter --snr 0.8 --fwhm 1.2 --background 1.5
```

Jeśli użyto `--filter`, ale nie podano żadnego progu, program zastosuje domyślne filtry jakości wzorca, tła i jasności gwiazd.

### 4.2. Progi bezwzględne

Opcje z sufiksem `-absolute` interpretują podaną liczbę jako rzeczywisty próg, a nie mnożnik mediany:

```bash
lucky-star ./sesja --filter --snr-absolute 10
```

Dostępne progi bezwzględne:

```text
--quality-absolute VALUE
--fwhm-absolute VALUE
--quality-star-pattern-absolute VALUE
--background-absolute VALUE
--star-brightness-absolute VALUE
--snr-absolute VALUE
```

Dla jednej metryki nie wolno jednocześnie podawać wersji względnej i bezwzględnej, np. `--snr 0.8 --snr-absolute 10`.

### 4.3. Folder `removed_*`

Odrzucone pliki są przenoszone do folderu, którego nazwa opisuje filtr, np.:

```text
removed_snr_10.000/
removed_fwhm_1.200/
removed_combined_1720000000/
```

Przy wielu filtrach używany jest folder `removed_combined_*`. Jeśli folder o danej nazwie już istnieje, program tworzy wariant z numerem, np. `removed_snr_10.000_2`.

## 5. Dzielenie najlepszych plików — `--divide`

`--divide` działa wyłącznie razem z `--filter` i przyjmuje dwa argumenty:

```text
--divide METRIC COUNT
```

Najpierw odrzucane pliki trafiają do `removed_*`. Następnie pozostałe pliki są sortowane według wskazanej metryki i przenoszone do kolejnych folderów zawierających maksymalnie `COUNT` plików.

### 5.1. Przykład z SNR

```bash
lucky-star ./sesja --filter --snr-absolute 5 --divide snr 1000
```

Jeżeli po filtrowaniu pozostało 3273 pliki, wynik będzie wyglądał następująco:

```text
1_snr/  # 1000 najlepszych plików
2_snr/  # kolejne 1000
3_snr/  # kolejne 1000
4_snr/  # ostatnie 273, najsłabsze z zachowanych
```

Dla `snr` większa wartość oznacza lepszy plik, więc sortowanie jest malejące.

### 5.2. Inne metryki

```bash
lucky-star ./sesja --filter --quality 0.8 --divide quality 500
lucky-star ./sesja --filter --divide fwhm 250
lucky-star ./sesja --filter --divide background 1000
lucky-star ./sesja --filter --divide star_brightness 1000
```

Dla `fwhm` i `background` mniejsza wartość jest lepsza, dlatego program sortuje je rosnąco. Dla pozostałych metryk sortowanie jest malejące.

Akceptowane nazwy metryk dla `--divide`:

```text
quality
fwhm
quality_star_pattern
background
star_brightness
snr
```

Można używać także myślników zamiast podkreśleń, np. `quality-star-pattern`.

Wszystkie pliki dzielone w ramach jednego uruchomienia muszą mieć dostępną wartość wybranej metryki. Program nie dzieli plików z brakującą wartością i zgłasza błąd.

## 6. Wzorzec gwiazd

Wzorzec pozwala porównywać te same gwiazdy między obrazami. Do utworzenia wzorca służy tryb interaktywny:

```bash
lucky-star --make-star-pattern ./sesja
```

Program utworzy obraz kandydatów, zwykle:

```text
star_pattern_candidates.jpg
```

Następnie należy wybrać numery gwiazd wskazane na obrazie:

```bash
lucky-star ./sesja --make-star-pattern ./sesja --star-pattern-numbers 1,4,7,9
```

Wzorzec zostanie zapisany jako:

```text
stars_pattern.json
```

Do analizy z użyciem wzorca należy podać jego ścieżkę:

```bash
lucky-star ./sesja \
  --star-pattern ./sesja/stars_pattern.json
```

Przykład filtrowania według metryk wzorca:

```bash
lucky-star ./sesja \
  --star-pattern ./sesja/stars_pattern.json \
  --filter \
  --snr 0.8 \
  --quality-star-pattern 0.9 \
  --divide snr 1000
```

## 7. Ograniczenie obszaru wyszukiwania gwiazd

Opcja `--crop FRACTION` ogranicza wyszukiwanie gwiazd do centralnej części obrazu. Przykład dla centralnych 30% szerokości i wysokości:

```bash
lucky-star ./sesja --crop 0.3
```

Można łączyć ją z filtrowaniem:

```bash
lucky-star ./sesja --crop 0.3 --filter --snr-absolute 8
```

## 8. Szybka kontrola najnowszych plików

`--check-count N` analizuje tylko `N` najnowszych plików. Tryb ten pomija wykresy i nie przenosi plików, dlatego służy do szybkiej kontroli:

```bash
lucky-star ./sesja --check-count 20
```

## 9. Zapisywanie obrazu z zaznaczonymi gwiazdami

Dla pojedynczego pliku można zapisać obraz z naniesionymi wykrytymi gwiazdami:

```bash
lucky-star ./sesja/klatka_0001.fits --save-stars
```

Program zapisuje pliki wynikowe obok analizowanego obrazu, w tym obraz JPG oraz informacje pomocnicze w formacie Markdown.

## 10. Plik konfiguracyjny

Program próbuje wczytać `config.json` z katalogu programu, a następnie z bieżącego katalogu. Jeśli plik nie istnieje lub nie można go odczytać, używane są wartości domyślne.

Przykładowy plik konfiguracyjny:

```json
{
  "gain_to_adu": {
    "5200": 27.5
  },
  "min_photons_to_detect_star": 10,
  "min_central_photons_to_detect_star": 20,
  "psf_size": 7,
  "min_photons_quality": 100.0,
  "rolling_avg_window": 10,
  "log_quality_window_t": 2.0,
  "star_pattern_position_tolerance_px": 5.0,
  "background_bias_adu": 0.0
}
```

Jeśli używasz gotowej konfiguracji projektu, najlepiej skopiować istniejący `config.json` i zmieniać tylko potrzebne wartości.

## 11. Typowy kompletny scenariusz

Poniższy przykład:

- używa wzorca gwiazd,
- analizuje centralną część obrazu,
- odrzuca obrazy zbyt słabe według kilku kryteriów,
- sortuje zachowane obrazy według SNR,
- tworzy grupy po 1000 plików.

```bash
lucky-star ./sesja \
  --star-pattern ./sesja/stars_pattern.json \
  --crop 0.5 \
  --filter \
  --snr 0.8 \
  --fwhm 1.2 \
  --background 1.5 \
  --divide snr 1000
```

Po zakończeniu:

```text
sesja/
├── removed_combined_*/
├── 1_snr/
├── 2_snr/
├── 3_snr/
├── ...
├── metrics_snr.png
├── metrics_fwhm.png
└── quality_map.*
```

Folder `1_snr` zawiera najlepsze zachowane obrazy, a kolejne foldery — następne grupy w kolejności jakości.

## 12. Ważne uwagi

- Operacje filtrowania i dzielenia przenoszą pliki, a nie tworzą ich kopii.
- Przed użyciem `--filter` warto wykonać najpierw analizę bez tej opcji.
- `--divide` wymaga `--filter`.
- `--divide` nie działa w trybie `--check-count`.
- Pliki niebędące FITS nie są analizowane.
- Metryki wzorca wymagają opcji `--star-pattern` oraz poprawnego dopasowania wzorca.
- Dla powtarzalności i bezpieczeństwa warto pracować na kopii sesji, szczególnie przy pierwszym użyciu filtrów.

## 13. Testy projektu

Uruchomienie testów jednostkowych:

```bash
cargo test
```

Kontrola kompilacji wersji produkcyjnej:

```bash
cargo build --release
```
