## 2024-05-30 - Unbounded file read and false positive vulnerability
**Vulnerability:** Unbounded file length read in TinyXML causing potential `bad_alloc` or buffer overflow when loading maliciously crafted 2GB+ files. The scanner flagged an `fgets` block that was already commented out.
**Learning:** Always check the memory allocations based on external input lengths, and remove unused or commented-out vulnerable patterns (like `fgets`) that may trigger false positive warnings in security tools.
**Prevention:** Impose strict, reasonable upper bounds (e.g., 500MB) on file reading operations where the whole file is slurped into memory.
