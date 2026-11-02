rem SPDX-FileCopyrightText: 2026 RFI
rem
rem SPDX-License-Identifier: Apache-2.0

set FILE=main
del %FILE%.py.lprof
uv run kernprof %FILE%.py -- --interactive -c config.toml
uv run python -m line_profiler -rmzt %FILE%.py.lprof > %FILE%.py.lines
echo Done. Read results in %FILE%.py.lines
