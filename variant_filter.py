# Process CSV into range objects
# AMP-LIST from Amp.csv file -> process_to_drop() -> assign_amps() -> count_amps() on master and drop list -> compare_amps()
import csv

class AmpliconRange(object):
    def __init__(self, start, stop, amp_id):
        self.range_ = (int(start), int(stop))
        self.amp_id = amp_id
        self.var_count = []

    def in_range(self, base_number):
        """
        Checks if base_number is in self.range_ from (start) to (stop) --> stop inclusive!
        """
        return (base_number in range(self.range_[0], self.range_[1]) )

    def add_variant(self, pos, chrom):
        self.var_count.append( (pos, chrom) )

    @staticmethod
    def query_amp(amp_id, amp_set):
        ret = None

        for item in amp_set:
            if amp_id == item.amp_id:
                ret = item
                break

        return ret

    def __len__(self):
        return len(self.var_count)

    def __repr__(self):
        return str(self.amp_id)


def read_amp_file(amp_file):
    range_objs = []
    with open(amp_file, 'rb') as amp_rdr:
        rdr = csv.DictReader(amp_rdr)
        for row in rdr:
            range_objs.append(AmpliconRange(
                amp_id=row['Amplicon_ID'],
                start=row['Amplicon_Start'],
                stop=row['Amplicon_Stop']
            ))

    return range_objs


def count_amps(amp_set, vcf_file):
    amp_dict = dict( (amp, 0) for amp in amp_set)

    for index in vcf_file.index:
        current_var = vcf_file.loc[index]
        current_amp = AmpliconRange.query_amp(
            amp_id=current_var['AMP'],
            amp_set=amp_set
        )
        testvar= amp_dict[current_amp]

        amp_dict[current_amp] += 1

    amp_dict = { amp: count for amp, count in amp_dict.items() if count }

    return amp_dict


def assign_amps(df, amps):
    """
    I hate pandas. i think its a great piece of whatever which is why everyone must be using it but who writes syntax like this ... who makes a data library so unintuitive and difficult to iterate through and modify values in. who makesit so hard to index specific columns in rows? i hate it. i want everyone to know i hate it. let it be known. Why does the documentation not have any examples at all. i'm stupid i cant figure this out. thats the whole reason that im at the documentation page right now. is everyone else really smart or just really quiet about how frustrating this library is.
    HALF the indexing and slicing methods dont work HALF THE time. every five seconds the type of the data youre working on will change into something else. last time a series type frame turned into a long long long short. thats really what happened. no slice methods worked on this at all. i.__dict__ returned an error. it devolved into a primative type
    Contact: slin63@illinois.edu
    """
    df = df.assign(amp=None)
    for index in df.index:
        where = tuple(df.loc[index][['CHR', 'POS']])
        amp_id = amp_sort(where[0], where[1], amps)
        df.loc[index, 'AMP'] = amp_id

    return df


def amp_sort(chr, pos, amps):
    """get hype for this N^2 run time boys"""
    for amp in amps:
        if amp.in_range(pos):
            ret = amp.amp_id
            break

    return ret


def compare_amps(master_dict, drop_dict, drop_list, master_list, cutoff_val):
    """
    Adds other variants to drop_list if >x% of variants on amplicon are on the drop list.
    Returns: List of amplicons to completely fail
    """
    to_drop = []
    # List comprehensions for optimal gains
    for amp in drop_dict:
        if ( float(drop_dict[amp]) / master_dict[amp] ) <= cutoff_val:
            to_drop.append(amp)

    return to_drop


def append_drop_list(drop_list, master_list, to_drop):
    """
    drop_list   :: DF with variants on drop_list
    master_list :: DF with all variants
    to_drop     :: list of amplicons that we should be dropping
    """
    for amp_site in to_drop:
        master_slice = master_list[master_list['AMP'] == amp_site.amp_id][['#CHROM', 'POS', 'AMP']]
        master_slice.columns = ['CHR', 'POS', 'AMP']
        drop_list = drop_list.append(master_slice)

    return drop_list





def process_to_drop(drop_list, master_list, amps, cutoff_val):
    """df = pandas.dataframe type object for the poor souls trying to understand this in the future"""
    to_assign = [drop_list, master_list]

    # Assign amplicons to each of the dataframes
    for i in range(len(to_assign)):
        to_assign[i] = assign_amps(to_assign[i], amps)

    # Count how many variants in the master_list are part of each amplicon
    amp_dict_master = count_amps(amps, to_assign[1])
    amp_dict_drop_l = count_amps(amps, to_assign[0])

    to_drop = compare_amps(
        master_dict=amp_dict_master,
        drop_dict=amp_dict_drop_l,
        drop_list=to_assign[0],
        master_list=to_assign[1],
        cutoff_val=cutoff_val
    )

    drop_list_appended = append_drop_list(
        drop_list=to_assign[0],
        master_list=to_assign[1],
        to_drop=to_drop
    )


    ## TODO: Filtering stuff tiffany asked for.

    return drop_list_appended, to_drop


# Debugging
if __name__ == '__main__':
    amp = 'ampliconregions.csv'
    read_amp_file(amp)