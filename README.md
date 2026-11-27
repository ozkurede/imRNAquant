# imRNAquant v1.0.0 by Ozkurede (2026)

imRNAquant quantifies reads containing uncleaved mitochondrial RNA junctions in
RNA-seq data. It searches for exact matches to 30-nucleotide uncleaved junction
identifiers (UJIs) and reports both the raw counts and counts normalized to
nearby mitochondrial transcript segments.

## Quick start

Use Python 3. The bundled example was tested on Linux with Python 3.12.14.
imRNAquant uses the Python standard library. There are no third-party packages to
install.

Open a terminal in the repository root, where `imRNAquant.py` and `CONFIG.py`
are located, and run:

```sh
python -c "from pathlib import Path; Path('output').mkdir(exist_ok=True)"
python imRNAquant.py
```

Use `python3` if that is the command for Python 3 on your system.
Press Enter after reviewing the settings displayed by the program.

The default `CONFIG.py` analyzes the bundled toy FASTQ files in
`example/fastq/` using the mouse reference. The run creates a timestamped log
and CSV report in `output/`. In the absolute-count table, the `M–Nd2`
junction should have 20 matches in Sample5 and Sample6 and 10 in each of the
other six samples.

## Analyze your own samples

Edit the root-level `CONFIG.py`. A copy of the toy configuration is available
at `example/CONFIG.py`.

1. Set `readFileFolder` to the directory containing your `.fastq.gz` files.
   Relative paths are resolved from the directory where you run the script.
2. Set `mtDseq_file`, `mtDmap_file`, and `mtD_UJI_file` to a matching set
   of FASTA, annotation JSON, and UJI JSON files from `resources/`.
3. Add one entry per sample to `samples_`. Give each sample a unique ID and
   assign it to the `Control` or `Test` group. At least three control samples are required.
4. Set `MC_max` to the number of worker processes appropriate for your
   machine. The default `read_process_limit` is 50 million reads per FASTQ
   file.

Use the orientation of the read sequences relative to the RNA transcript to
choose the FASTQ fields. Check the library preparation protocol to determine
the orientation. Filename suffixes such as `_R1` and `_R2` identify read
numbers.

| Input | `reverse` | `forward` |
| --- | --- | --- |
| Paired-end, R1 reverse and R2 forward | `"sample_R1.fastq.gz"` | `"sample_R2.fastq.gz"` |
| Single-end, reverse-oriented | `"sample_R1.fastq.gz"` | `None` |
| Single-end, forward-oriented | `None` | `"sample_R1.fastq.gz"` |

For example, a forward-oriented single-end control sample can be entered as:

```python
{
    "id": "Sample1",
    "group": "Control",
    "reverse": None,
    "forward": "sample_R1.fastq.gz",
},
```

## Output

- `*_log.txt` records configuration settings and run progress.
- `*_report.csv` contains absolute read counts, junction counts per 1,000
  reads matching the corresponding reference segment, and values relative
  to the control mean. `Ref0` marks a zero reference-segment count, and
  `!LowCt` marks a result affected by low counts in the controls. `Unavail` 
  indicates that fold changes were unavailable because at least one control
  had zero reads matching the corresponding reference segment.

  The final two sections of the `*_report.csv` file contain the main imRNA results.
    - **imRNAs retaining 5′ tRNA extensions** appear under the heading
      `Reporting relative count for jType tR-mR normalized to refPost_head`.
    - **imRNAs retaining 3′ tRNA extensions** appear under the heading
      `Reporting relative count for jType mR-tR normalized to refPre_tail`.

    Each section contains two tables placed side by side, separated by an empty
    column. The left table reports junction counts per 1,000 reads matching the
    corresponding reference segment. The right table expresses each sample’s
    relative count as a multiple of the control mean, so a value of 1 equals
    the control mean.

## Bundled resources

The human and mouse mitochondrial reference sequences are identified in
their FASTA headers by RefSeq accessions `NC_012920.1` and `NC_005089.1`,
respectively. The `mouse_mtDNA_6NCrl_*` files provide a separate C57BL/6NCrl
reference set. The FASTQ files in `example/fastq/` are a toy dataset for
checking that the workflow runs and produces the expected counts.

## License

The imRNAquant source code is available under the MIT License. See
[LICENCE](LICENCE).
