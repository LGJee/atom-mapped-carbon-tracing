"""Disk-backed propagation states with bounded, ordered substrate reads.

Only the parent allocates IDs and merges states, at the existing layer/round
boundaries. A product discovered during a round is never a substrate in that
same round. SQLite files are temporary working data, not resumable checkpoints.
"""
from collections import OrderedDict
from functools import lru_cache
import json

PAGE_SIZE = 256
CACHE_PAGES = 8
CATALOG_CACHE = 4096


class DiskCatalog:
    def __init__(self, db):
        self.db, self.count = db, 0
        self.cache = OrderedDict()

    def remember(self, key, number):
        self.cache[key] = number
        self.cache.move_to_end(key)
        if len(self.cache) > CATALOG_CACHE:
            self.cache.popitem(last=False)

    def lookup(self, name, marked):
        key = (name, marked)
        if key in self.cache:
            self.cache.move_to_end(key)
            return self.cache[key]
        row = self.db.execute('SELECT id FROM catalog WHERE name=? AND marked=?', key).fetchone()
        if row:
            self.remember(key, row[0])
            return row[0]
        return None

    def __contains__(self, key):
        return self.lookup(*key) is not None

    def __getitem__(self, key):
        value = self.lookup(*key)
        if value is None:
            raise KeyError(key)
        return value

    def __setitem__(self, key, number):
        self.db.execute('INSERT INTO catalog VALUES (?, ?, ?)', (number,) + key)
        self.count += 1
        self.remember(key, number)

    def __len__(self):
        return self.count


class DiskSteps:
    def __init__(self, db, states):
        self.db, self.states, self.revision = db, states, 0

    def get(self, key, default=None):
        row = self.db.execute('SELECT step FROM states WHERE name=? AND id=?', key).fetchone()
        return row[0] if row and row[0] is not None else default

    def __getitem__(self, key):
        value = self.get(key)
        if value is None:
            raise KeyError(key)
        return value

    def __setitem__(self, key, step):
        previous = self.get(key)
        if previous != step:
            self.db.execute('INSERT INTO states(name,id,step) VALUES (?,?,?) '
                'ON CONFLICT(name,id) DO UPDATE SET step=excluded.step', key + (step,))
            self.revision += 1
            self.states.page.cache_clear()

    def items(self):
        for name, number, step in self.db.execute('SELECT name,id,step FROM states WHERE step IS NOT NULL'):
            yield (name, number), step


class MetaboliteStates:
    def __init__(self, owner, name):
        self.owner, self.name = owner, name

    def __len__(self):
        return self.owner.counts.get(self.name, 0)

    def __contains__(self, number):
        return self.owner.db.execute('SELECT 1 FROM states WHERE name=? AND id=? AND assignment IS NOT NULL',
                                     (self.name, number)).fetchone() is not None

    def __setitem__(self, number, assignment):
        exists = number in self
        self.owner.db.execute('INSERT INTO states(name,id,assignment) VALUES (?,?,?) '
            'ON CONFLICT(name,id) DO UPDATE SET assignment=excluded.assignment',
            (self.name, number, json.dumps(assignment, separators=(',', ':'))))
        if not exists:
            self.owner.counts[self.name] = len(self) + 1
        self.owner.page.cache_clear()

    def items(self):
        after = -1
        for _ in range(0, len(self), PAGE_SIZE):
            page = self.owner.page(self.name, after)
            for number, assignment, _ in page:
                yield number, assignment
            after = page[-1][0]

    def values(self):
        for _, assignment in self.items():
            yield assignment


class DiskStates:
    def __init__(self, db, memory_guard):
        self.db, self.guard, self.counts = db, memory_guard, {}
        self.page = lru_cache(maxsize=CACHE_PAGES)(self._page)

    def _page(self, name, after):
        if self.guard:
            size = self.db.execute('SELECT COALESCE(SUM(length(assignment)),0) FROM '
                '(SELECT assignment FROM states WHERE name=? AND assignment IS NOT NULL AND id>? '
                'ORDER BY id LIMIT ?)', (name, after, PAGE_SIZE)).fetchone()[0]
            self.guard.check('reading a propagation-state page', force=True, extra_bytes=size * 16)
        rows = self.db.execute('SELECT id,assignment,step FROM states '
            'WHERE name=? AND assignment IS NOT NULL AND id>? ORDER BY id LIMIT ?',
            (name, after, PAGE_SIZE)).fetchall()
        return tuple((number, tuple(tuple(a) for a in json.loads(assignment)), step)
                     for number, assignment, step in rows)

    def __getitem__(self, name):
        return MetaboliteStates(self, name)

    def get(self, name, default=None):
        return self[name] if self.counts.get(name, 0) else default

    def items(self):
        for name in self.counts:
            yield name, self[name]

    def choices(self, groups, required_labels):
        return [StateChoices(self, name, transfers, required_labels) for _, name, transfers in groups]


class StateChoices:
    def __init__(self, states, name, transfers, required):
        self.states, self.name, self.transfers, self.required = states, name, transfers, required
        self.count = states.counts.get(name, 0)

    def __len__(self):
        return self.count or 1

    def __iter__(self):
        if not self.count:
            yield None
            return
        after = -1
        for _ in range(0, self.count, PAGE_SIZE):
            page = self.states.page(self.name, after)
            for number, assignment, first in page:
                labels = {(part, atom): label for part, atom, label in assignment}
                yield (number, tuple((mapped, labels[(part, atom)])
                       for part, atom, mapped in self.transfers if (part, atom) in labels), first,
                       self.required is None or any(label in self.required for _, _, label in assignment))
            after = page[-1][0]


def ordered_product(choices):
    """Match itertools.product order without materializing each substrate pool."""
    if not choices:
        yield ()
        return
    iterators, selected = [iter(choices[0])], []
    while iterators:
        try:
            value = next(iterators[-1])
        except StopIteration:
            iterators.pop()
            if selected:
                selected.pop()
            continue
        if len(iterators) == len(choices):
            yield tuple(selected + [value])
        else:
            selected.append(value)
            iterators.append(iter(choices[len(iterators)]))


def create_state_tables(db, guard):
    db.execute('CREATE TABLE catalog(id INTEGER PRIMARY KEY, name TEXT, marked TEXT, UNIQUE(name,marked))')
    db.execute('CREATE TABLE states(name TEXT, id INTEGER, assignment TEXT, step INTEGER, PRIMARY KEY(name,id)) WITHOUT ROWID')
    states = DiskStates(db, guard)
    return states, DiskCatalog(db), DiskSteps(db, states)
