print(f"Loading {__file__!r} ...")

def my_export(str_list, interval=1,det = 'merlin1', mon = 'sclr1_ch4'):
    sid_list = get_sid_list(str_list = str_list,interval=interval)
    num = np.size(sid_list)
    for sid in sid_list:
        #sid, df = _load_scan(sid, fill_events=False)
        sid = int(sid)
        h = db[sid]
        #sid = h.start['scan_id']
        df = h.table()
        mots = h.start['motors']

        dir = os.path.join('/data/home/home/hyan/export','scan_{:06d}'.format((sid//10000)*10000))
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
            path = os.path.join(dir, 'scan_{}_{}.h5'.format(sid,det))
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
            path = os.path.join(dir, 'scan_{}_{}.h5'.format(sid,det))
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

def my_export_diffraction(str_list, interval, name_list,det = 'merlin1'):
    sid_list = get_sid_list(str_list = str_list,interval=interval)
    num = np.size(sid_list)
    
    for sid in sid_list:
        sid = int(sid)
        hdr = db[sid]
        df = hdr.table()
        #sid = hdr.start['scan_id']
        dir = os.path.join('/data/home/home/hyan/export','scan_{:06d}'.format((sid//10000)*10000))
        if os.path.exists(dir) == False:
            print('{} does not exist.'.format(dir))
            mycmd = ''.join(['mkdir',' ',dir])
            os.system(mycmd)
            if os.path.exists(dir):
                print('{} created successfully'.format(dir))
            else:
                print('Can''t create {}. Quit exporting '.format(dir))
                return
        path = os.path.join(dir, 'scan_{}_{}.txt'.format(sid,det))
        df.to_csv(path, float_format='%1.5e', sep='\t', columns=name_list)
        print('Scan {}. Saving to {}'.format(sid, path))
        images = list(db[sid].data(det))
        images = np.squeeze(images)
        path = os.path.join(dir, 'scan_{}_{}.h5'.format(sid,det))
        f = h5py.File(path, 'w')
        dset = f.create_dataset('/entry/instrument/detector/data', data=images)
        f.close()
        print('Scan {}. Saving to {}'.format(sid, path))



def my_export_angle(str_list,interval,rot_motor):
    j = 0
    sid_list = get_sid_list(str_list = str_list,interval=interval)
    num = np.size(sid_list)
    #num = round((sid_end-sid_start)/interval)
    angle = zeros((num,2))
    for sid in sid_list:
        sid = int(sid)
        dir = os.path.join('/data/home/home/hyan/export','scan_{:06d}'.format((sid//10000)*10000))
        if os.path.exists(dir) == False:
            print('{} does not exist.'.format(dir))
            mycmd = ''.join(['mkdir',' ',dir])
            os.system(mycmd)
            if os.path.exists(dir):
                print('{} created successfully'.format(dir))
            else:
                print('Can''t create {}. Quit exporting '.format(dir))
                return
        angle[j,0] = sid
        angle[j,1] = check_baseline(sid,rot_motor) 
        j = j+1 
    path = os.path.join(dir, 'scan_{}_{}_angle.txt'.format(sid_list[0],sid_list[-1]))
    np.savetxt(path, angle, fmt ='%d %10.5f')
    print('Angle list for scan {} to {} saved to {}'.format(sid_list[0],sid_list[1], path))    

def get_sid_list(str_list, interval):
    num_elem = np.size(str_list)
    for i in range(num_elem):
        str_elem = str_list[i].split('-')
        if i == 0:
            if np.size(str_elem) == 1:
                tmp = int(str_elem[0])
            else:
                tmp = np.arange(int(str_elem[0]),int(str_elem[1])+1,interval)
            sid_list = np.reshape(tmp,(-1,))
        else:
            if np.size(str_elem) == 1:
                tmp = int(str_elem[0])
            else:
                tmp = np.arange(int(str_elem[0]),int(str_elem[1])+1,interval)
            tmp = np.reshape(tmp,(-1,))
            sid_list = np.concatenate((sid_list,tmp))
    return sid_list
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
