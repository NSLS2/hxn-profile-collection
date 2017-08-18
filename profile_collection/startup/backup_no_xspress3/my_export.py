def my_export(sid,num=1):
    for i in range(num):
        sid, df = _load_scan(sid, fill_events=False)
        path = os.path.join('/home/hyan/export/', 'scan_{}.txt'.format(sid))
        print('Scan {}. Saving to {}'.format(sid, path))
        #non_objects = [name for name, col in df.iteritems() if col.dtype.name not in ('object', )]
        #dump all data
        non_objects = [name for name, col in df.iteritems()]
        df.to_csv(path, float_format='%1.5e', sep='\t',
                  columns=sorted(non_objects))

        path = os.path.join('/home/hyan/export/', 'scan_{}_scaler.txt'.format(sid))
        #np.savetxt(path, (df['sclr1_ch3'], df['p_ssx'], df['p_ssy']), fmt='%1.5e')
        np.savetxt(path, (df['sclr1_ch4'], df['zpssx'], df['zpssy']), fmt='%1.5e')
        path = os.path.join('/home/hyan/export/', 'scan_{}.h5'.format(sid))
        filename = get_all_filenames(sid,'merlin1')
        for fn in filename:
            break
        mycmd = ''.join(['scp', ' ', fn, ' ', path])
        os.system(mycmd)

        sid = sid + 1

