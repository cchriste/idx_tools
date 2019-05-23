```
range - incrementally computes range of fields in a .idx volume
Copyright (C) 2019  Cameron Christensen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
```

# range

* incrementally computes range of fields in a .idx volume


## Usage

range [params] <dataset>

### dataset

* location of dataset, either local or remote

### params

-a,--all            compute range of all fields in volume

-f,--field <name>   field for which range should be computed

-l,--list           list all fields and their current ranges

-o,--output <name>  output idx path, by default same as <dataset>

--min/--max         ignore any values smaller/larger than these
                    some datasets use placeholders for default/nil values

-v,--verbose        show more output of the task in progress

