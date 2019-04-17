

def my_export(sid,num=1, interval=1,det = 'merlin1', mon = 'sclr1_ch4'):
    for i in range(num):
        #sid, df = _load_scan(sid, fill_events=False)
        h = db[sid]
        sid = h.start['scan_id']
        df = h.table()
        mots = h.start['motors']

        dir = os.path.join('/data/home/hyan/export','scan_{:06d}'.format((sid//10000)*10000))
        if os.path.exists(dir) == False:
            print('{} does not exist.'.format(dir))
            mycmd = ''.join(['mkdir',' ',dir])
            os.system(mycmd)
            if os.path.exists(dir):
                print('{} created successfully'.format(dir))
            else:
                print('Can''t create {}. Quit exporting '.format(dir))
                return
        path = os.path.join(dir, 'scan_{}.txt'.format(sid))
        print('Scan {}. Saving to {}'.format(sid, path))
        #non_objects = [name for name, col in df.iteritems() if col.dtype.name not in ('object', )]
        #dump all data
        non_objects = [name for name, col in df.iteritems()]
        df.to_csv(path, float_format='%1.5e', sep='\t',
                  columns=sorted(non_objects))
        path = os.path.join(dir,'scan_{}_scaler.txt'.format(sid))
        #np.savetxt(path, (df['sclr1_ch3'], df['p_ssx'], df['p_ssy']), fmt='%1.5e')
        np.savetxt(path, (df[mon], df[mots[0]], df[mots[1]]), fmt='%1.5e')

        filename = get_path(sid,det)
        num_subscan = len(filename)
        #num_subscan = 2
        if num_subscan == 1:
            for fn in filename:
                break
            path = os.path.join(dir, 'scan_{}.h5'.format(sid))
            mycmd = ''.join(['scp', ' ', fn, ' ', path])
            os.system(mycmd)
            #num_subscan=-1
        else:
            #h = db[sid]
            #df = db.get_table(h,fill=False)
            images = list(db[sid].data(det))
            imgs = np.squeeze(images)
            '''
            num_frame, tmp = np.shape(df)
            for i in range(num_frame):
                if np.mod(i,1000) == 0:
                    print('load frame ', i)
                image = np.squeeze(filestore.api.get_data(df['merlin1'][i+1]))
                if i == 0:
                    nx,ny = np.shape(image)
                    images = np.zeros((num_frame,nx,ny))
                images[i,:,:] = image
            '''
            path = os.path.join(dir, 'scan_{}.h5'.format(sid))
            f = h5py.File(path, 'w')
            dset = f.create_dataset('/entry/instrument/detector/data', data=imgs)
            f.close()

        print('Scan {}. Saving to {}'.format(sid, path))


        '''
            j = 1
            for fn in filename:
                path = os.path.join('/data/home/hyan/export/', 'scan_{}_{}.h5'.format(sid, j))
                mycmd = ''.join(['scp', ' ', fn, ' ', path])
                os.system(mycmd)
                j = j + 1
        '''
        sid = sid + interval

def my_export_1d(sid_start, sid_end, name_list, interval = 1, det = 'merlin1'):
    for i in range (sid_start, sid_end+1,interval):
        hdr = db[i]
        df = hdr.table()
        sid = hdr.start['scan_id']
        dir = os.path.join('/data/home/hyan/export','scan_{:06d}'.format((sid//10000)*10000))
        if os.path.exists(dir) == False:
            print('{} does not exist.'.format(dir))
            mycmd = ''.join(['mkdir',' ',dir])
            os.system(mycmd)
            if os.path.exists(dir):
                print('{} created successfully'.format(dir))
            else:
                print('Can''t create {}. Quit exporting '.format(dir))
                return
        path = os.path.join(dir, 'scan_{}.txt'.format(sid))
        df.to_csv(path, float_format='%1.5e', sep='\t', columns=name_list)
        print('Scan {}. Saving to {}'.format(sid, path))
        images = list(db[sid].data(det))
        images = np.squeeze(images)
        path = os.path.join(dir, 'scan_{}.h5'.format(sid))
        f = h5py.File(path, 'w')
        dset = f.create_dataset('/entry/instrument/detector/data', data=images)
        f.close()
        print('Scan {}. Saving to {}'.format(sid, path))





'''
def get_all_filenames(scan_id, key='merlin1'):
    #scan_id, df = _load_scan(scan_id, fill_events=False)
    h = db[scan_id]
    df = db.get_table(h)
    from filestore.path_only_handlers import (AreaDetectorTiffPathOnlyHandler,
                                              RawHandler)
    handlers = {'AD_TIFF': AreaDetectorTiffPathOnlyHandler,
                'XSP3': RawHandler,
                'AD_HDF5': RawHandler,
                'TPX_HDF5': RawHandler,
                }
    # this is easy to change under new databroker
    # hdr = db[scan_id], hdr.db.reg.
    filenames = [filestore.api.retrieve(uid, handlers)[0]
                 for uid in list(df[key])]

    if len(set(filenames)) != len(filenames):
        return set(filenames)
    return filenames
'''
