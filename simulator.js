/** @typedef {{x:number, y:number, z:number}} Vector3 */

/**
 * @typedef {Object} PointData - 頂点１つのデータ
 * @property {number} index
 * @property {Vector3[]} triangle
 * @property {Vector3} centerOfGravity
 * @property {Vector3} normalVector
 * @property {{point:Vector3, normalVector:Vector3}[]} insideJudge
 */
/**
 * オブジェクトの頂点データ（衝突判定で使用）
 * @type {{normalVectors:{data:PointData[], sort:{x:PointData[], y:PointData[], z:PointData[]}}}}
 */
let data = {};

/**
 * 加速度を計算するときに使う項
 * @type {{pressureTerm:Vector3, viscosityTerm:Vector3}[]}
 */
const terms = [];

const h = 0.3/* 0.012 */; //影響半径
const particleMass = 0.0002; //粒子の質量

const g = createVector3(0, -9.8, 0);  // 重力加速度
const pressureStiffness = 200; //圧力係数
const restDensity = 1000; //静止密度
const viscosity = 0.000001;  // 粘性係数
const attenuationCoefficient = -5;  // ダンパ係数
const springConstant = -5;  // ばね係数

const densityCoef = particleMass * 315 / (64 * Math.PI * Math.pow(h,9)); //密度計算で使う

const pressureCoef = particleMass * 45 / (Math.PI * Math.pow(h,6)); //圧力項計算で使う
const viscosityCoef = viscosity * particleMass * 45 / (Math.PI * Math.pow(h,6)); //粘性項計算で使う

//壁は立方体
const wallWidth = 10; 
const thickness = 3/* 0.5, 0.1 */; //壁となる粒子の厚み（個数）
const interval = h / thickness;


const shaft = ["x","y","z"];

//壁となる、速度がゼロで固定の粒子群
function makeWall(particles){
    let newPosition = createVector3(0,0,0);
    for(let x = 0; x < 11 ; x = x + 10){
        newPosition[shaft[2]] = x;
        for(let i=0 ; i< wallWidth / interval; i++){//高さ
            for(let j = 0; j < wallWidth / interval; j++){//横幅
                newPosition[shaft[1]] = (interval / 2) + i * interval;
                newPosition[shaft[0]] = (interval / 2) + j * interval;
                const newParticle = createParticle(newPosition, true);
                particles.push(newParticle);
            }

        }  
    }
    for(let x = 0; x < 11 ; x = x + 10){
        newPosition[shaft[0]] = x;
        for(let i=0 ; i< wallWidth / interval; i++){//高さ
            for(let j = 0; j < wallWidth / interval; j++){//横幅
                newPosition[shaft[1]] = (interval / 2) + i * interval;
                newPosition[shaft[2]] = (interval / 2) + j * interval;
                const newParticle = createParticle(newPosition, true);
                particles.push(newParticle);
            }

        }  
    }
    newPosition[shaft[1]] = 0;
        for(let i=0 ; i< wallWidth / interval; i++){//x方向
            for(let j = 0; j < wallWidth / interval; j++){//z方向
                newPosition[shaft[0]] = (interval / 2) + i * interval;
                newPosition[shaft[2]] = (interval / 2) + j * interval;
                const newParticle = createParticle(newPosition, true);
                particles.push(newParticle);
            }

        }
    //console.log(particles)
}
function addParticles(particles) {
    for (let i = 2; i < 5000; i += 2) {
        particles.push(createParticle(createVector3(Math.random() + 0.5, i, 1)));
    }
}



/**
 * 粒子の密度計算
 * @param {{position:Vector3, velocity:Vector3, force:Vector3, density:number, pressure:number}[]} particles - 粒子のリスト
 */
function calcDensity(particles) {
    const h2 = h*h; //事前にhの二乗を計算しておく
    for (let i = 0; i < particles.length; i++) { //一つづつ粒子の密度を計算

        reportProgress(reportHeader + "密度を計算中..." + i + "/" + particles.length);

        let nowParticle = particles[i]; //今回計算する粒子
        let sum = 0; //足し合わせる変数
        for (let j = 0; j < particles.length; j++) { //他の粒子全てについて
            if(i == j){continue;} //自分自身だったらスキップ
            let nearParticle = particles[j];
            
            let diff = subVector3(nearParticle.position, nowParticle.position); //粒子距離
            const r2 = dotVector3(diff, diff); //粒子距離の２乗

            //粒子距離がhより小さい場合だけ計算する
            if ( r2 < h2 ) {
                const c = h2 - r2;
                sum += Math.pow(c,3); //(h2-r2)の３乗
            }
        }

        nowParticle.density = sum * densityCoef; //密度が求まった
    }
}

/**
 * 粒子の圧力計算
 * @param {{position:Vector3, velocity:Vector3, force:Vector3, density:number, pressure:number}[]} particles - 粒子のリスト
 */
function calcPressure(particles) {
    for(let i = 0; i < particles.length; i++) { //一つづつ粒子の圧力を計算

        reportProgress(reportHeader + "圧力を計算中..." + i + "/" + particles.length);

        particles[i].pressure = pressureStiffness * ( particles[i].density - restDensity );
    }
}

/**
 * 粒子の圧力項計算
 * @param {{position:Vector3, velocity:Vector3, force:Vector3, density:number, pressure:number}[]} particles - 粒子のリスト
 */
function calcPressureTerm(particles) {
    const h2 = h*h; //事前にhの二乗を計算しておく
    for (let i = 0; i < particles.length; i++) { //一つづつ粒子の密度を計算

        reportProgress(reportHeader + "圧力項を計算中..." + i + "/" + particles.length);

        let nowParticle = particles[i]; //今回計算する粒子
        if (nowParticle.is_wall) continue;  // 速度の計算や位置の更新をしない粒子の場合スキップ
        let sum = createVector3(0, 0, 0); //足し合わせる変数
        for (let j = 0; j < particles.length; j++) { //他の粒子全てについて
            if(i == j){continue;} //自分自身だったらスキップ
            let nearParticle = particles[j];
            
            let diff = subVector3(nearParticle.position, nowParticle.position); //粒子距離
            const r2 = dotVector3(diff, diff); //粒子距離の２乗

            //粒子距離がhより小さい場合だけ計算する
            if ( r2 < h2 ) {
                const r = Math.sqrt(r2); //粒子距離
                const c = h - r;
                const n = ((nearParticle.pressure /*-*/+ nowParticle.pressure) / (2 * nearParticle.density)) * Math.pow(c,2) / r;
                sum = addVector3(sum, multiplyScalarVector3(diff, n));
            }
        }

        if (!terms[i]) terms[i] = {};
        terms[i].pressureTerm = multiplyScalarVector3(sum, (-1/*/nowParticle.pressure*/) * pressureCoef);  // 圧力項が求まった
    }
}

/**
 * 粒子の粘性項計算
 * @param {{position:Vector3, velocity:Vector3, force:Vector3, density:number, pressure:number}[]} particles - 粒子のリスト
 */
function calcViscosityTerm(particles) {
    const h2 = h*h; //事前にhの二乗を計算しておく
    for (let i = 0; i < particles.length; i++) { //一つづつ粒子の密度を計算

        reportProgress(reportHeader + "粘性項を計算中..." + i + "/" + particles.length);

        let nowParticle = particles[i]; //今回計算する粒子
        if (nowParticle.is_wall) continue;  // 速度の計算や位置の更新をしない粒子の場合スキップ
        let sum = createVector3(0, 0, 0); //足し合わせる変数
        for (let j = 0; j < particles.length; j++) { //他の粒子全てについて
            if(i == j){continue;} //自分自身だったらスキップ
            let nearParticle = particles[j];
            
            let diff = subVector3(nearParticle.position, nowParticle.position); //粒子距離
            const r2 = dotVector3(diff, diff); //粒子距離の２乗

            //粒子距離がhより小さい場合だけ計算する
            if ( r2 < h2 ) {
                const r = Math.sqrt(r2); //粒子距離
                const c = h - r;
                const n = c / nearParticle.density;
                sum = addVector3(sum, multiplyScalarVector3(subVector3(nearParticle.velocity, nowParticle.velocity), n));
            }
        }

        if (!terms[i]) terms[i] = {};
        terms[i].viscosityTerm = multiplyScalarVector3(sum, viscosityCoef);  // 粘性項が求まった
    }
}

/**
 * 範囲に合う面を抽出
 * @param {boolean[]} resultArray - 結果を出力する配列
 * @param {PointData[]} sortedPointArray - ソートされた面の配列
 * @param {number} min - 抽出する範囲の最小値
 * @param {number} max - 抽出する範囲の最大値
 * @param {"x" | "y" | "z"} axis - 抽出する軸
 */
function filterPoint(resultArray, sortedPointArray, min, max, axis) {
    let start = 0;
    let end = sortedPointArray.length - 1;
    let middle = start + Math.round((end - start) / 2);
    let count = 0;
    while (start > end && !(min <= sortedPointArray[middle].centerOfGravity[axis] && sortedPointArray[middle].centerOfGravity[axis] <= max)) {
        if (sortedPointArray[middle].centerOfGravity[axis] < min) {
            start = middle + 1;
        }
        if (max < sortedPointArray[middle].centerOfGravity[axis]) {
            end = middle - 1;
        }
        middle = start + Math.round((end - start) / 2);
        count++;
        if (count >= 10000) break;
    }

    let start_min = start;
    let end_min = middle;
    start = Math.floor((end_min - start_min) / 2);
    count = 0;
    while (start > 0 && !(sortedPointArray[start].centerOfGravity[axis] >= min && sortedPointArray[start-1].centerOfGravity[axis] < min)) {
        if (sortedPointArray[start].centerOfGravity[axis] < min) {
            start_min = start + 1;
        }
        if (sortedPointArray[start].centerOfGravity[axis] > min) {
            end_min = start - 1;
        }
        start = start_min + Math.floor((end_min - start_min) / 2);
        count++;
        if (count >= 10000) break;
    }

    let start_max = middle;
    let end_max = end;
    end = start_max + Math.ceil((end_max - start_max) / 2);
    count = 0;
    while (end < sortedPointArray.length - 1 && !(sortedPointArray[end].centerOfGravity[axis] <= max && sortedPointArray[end+1].centerOfGravity[axis] > max)) {
        if (sortedPointArray[end].centerOfGravity[axis] < max) {
            start_max = end + 1;
        }
        if (sortedPointArray[end].centerOfGravity[axis] > max) {
            end_max = end - 1;
        }
        end = start_max + Math.ceil((end_max - start_max) / 2);
        count++;
        if (count >= 10000) break;
    }

    for (let i = 0; i < start; i++) {
        resultArray[i] = undefined;
    }
    for (let i = start; i <= end; i++) {
        resultArray[i] = resultArray[i] && true;
    }
    for (let i = end; i < sortedPointArray.length; i++) {
        resultArray[i] = undefined;
    }
}
/**
 * 粒子の内在判定
 * @param {Vector3} position - 粒子の位置
 * @param {PointData} point - 判定する面
 */
function is_inside(position, point) {
    for (let i = 0; i < point.insideJudge.length; i++) {
        const dot = dotVector3(point.insideJudge[i].normalVector, subVector3(position, point.insideJudge[i].point));
        if (dot < 0) return false;
    }
    return true;
}
/**
 * 粒子の衝突項計算
 * @param {{position:Vector3, velocity:Vector3, force:Vector3, density:number, pressure:number}[]} particles - 粒子のリスト
 */
function calcColiderTerm(particles) {
    for (let i = 0; i < particles.length; i++) {

        reportProgress(reportHeader + "衝突項を計算中..." + i + "/" + particles.length);

        let nowParticle = particles[i]; //今回計算する粒子
        let term = createVector3();

        const filterArray = [];
        filterPoint(filterArray, data.normalVectors.sort.x, nowParticle.position.x - (-attenuationCoefficient), nowParticle.position.x + (-attenuationCoefficient), "x");
        filterPoint(filterArray, data.normalVectors.sort.y, nowParticle.position.y - (-attenuationCoefficient), nowParticle.position.y + (-attenuationCoefficient), "y");
        filterPoint(filterArray, data.normalVectors.sort.z, nowParticle.position.z - (-attenuationCoefficient), nowParticle.position.z + (-attenuationCoefficient), "z");
        for (let j = 0; j < data.normalVectors.data.length; j++) {
            if (!filterArray[j]) continue;
            const nowPoint = data.normalVectors.data[j];
            if (is_inside(nowParticle.position, nowPoint)) {
                let distance = Math.abs(nowPoint.normalVector.x * (nowParticle.position.x - nowPoint.centerOfGravity.x) + nowPoint.normalVector.y * (nowParticle.position.y - nowPoint.centerOfGravity.y) + nowPoint.normalVector.z * (nowParticle.position.z - nowPoint.centerOfGravity.z));
                let nowTerm = multiplyScalarVector3(nowPoint.normalVector, springConstant * distance + attenuationCoefficient * dotVector3(nowParticle.velocity, nowPoint.normalVector));
                term = addVector3(term, nowTerm);
            }
        }
        console.log(term);
        terms[i].coliderTerm = term;

        /*let distance = nowParticle.position.y;
        if (distance < -attenuationCoefficient) {
            terms[i].coliderTerm = multiplyScalarVector3(createVector3(0, 1, 0), springConstant * distance + attenuationCoefficient * dotVector3(nowParticle.velocity, createVector3(0, 1, 0)));
        } else {
            terms[i].coliderTerm = createVector3();
        }*/
    }
}


/**
 * 粒子のリスト
 * @type {{
 *      position:Vector3,
 *      velocity:Vector3,
 *      force:Vector3,
 *      acceleration:Vector3,
 *      density:number,
 *      pressure:number,
 *      is_wall:boolean
 * }[]}
 */
const _particles = [];


const deltaTime = 0.01;  // （秒）

let reportHeader = "";

function tick() {
    calcDensity(_particles);
    calcPressure(_particles);
    calcPressureTerm(_particles);
    calcViscosityTerm(_particles);
    calcColiderTerm(_particles);

    for (let i = 0; i < _particles.length; i++) {
        const nowParticle = _particles[i];
        if (nowParticle.is_wall) continue;
        const a = addVector3(addVector3(addVector3(terms[i].pressureTerm, terms[i].viscosityTerm), terms[i].coliderTerm), g);
        const v = addVector3(nowParticle.velocity, multiplyScalarVector3(addVector3(nowParticle.acceleration, a), 0.5 * deltaTime));
        const deltaPosition = addVector3(multiplyScalarVector3(v, deltaTime), multiplyScalarVector3(a, 0.5 * deltaTime * deltaTime));
        nowParticle.acceleration = a;
        nowParticle.velocity = v;
        nowParticle.position = addVector3(nowParticle.position, deltaPosition);
    }
}

/**
 * シミュレーションを開始する関数
 * @param {number} simulateSeconds - シミュレーションを行う秒数
 */
function start(simulateSeconds) {
    // 初めに実行する処理
    //makeWall(_particles);
    addParticles(_particles);

    self.postMessage({type:"result", content:_particles, time:0});

    for (let i = 0; i < simulateSeconds / deltaTime; i++) {
        const time = (i + 1) * deltaTime;
        reportHeader = time + "秒目/" + simulateSeconds + "秒 ";
        tick();
        self.postMessage({type:"result", content:_particles, time:time});
    }

    reportProgress("シミュレーション終了");
}
self.addEventListener("message", event => {
    data = event.data;
    for (let i = 0; i < data.normalVectors.data.length; i++) {
        data.normalVectors.data[i].insideJudge.push({
            point: addVector3(data.normalVectors.data[i].triangle[0], multiplyScalarVector3(data.normalVectors.data[i].normalVector, -attenuationCoefficient)),
            normalVector: multiplyScalarVector3(data.normalVectors.data[i].normalVector, -1)
        });
        data.normalVectors.data[i].insideJudge.push({
            point: subVector3(data.normalVectors.data[i].triangle[0], multiplyScalarVector3(data.normalVectors.data[i].normalVector, -attenuationCoefficient)),
            normalVector: data.normalVectors.data[i].normalVector
        });
    }
    start(100);
});



/* 関数定義部分 */

/**
 * 粒子のオブジェクトを作成する関数
 * @param {Vector3} position - 粒子の位置
 * @param {Vector3} velocity - 粒子の速度
 * @param {Vector3} force - 粒子にかかっている力
 * @param {Vector3} acceleration - 粒子の加速度
 * @param {number} density - 粒子の密度
 * @param {number} pressure - 粒子の圧力
 * @param {boolean} is_wall - 粒子が壁かどうか
 * @returns {{position:Vector3, velocity:Vector3, force:Vector3, acceleration:Vector3, density:number, pressure:number, is_wall:boolean}} - 粒子のオブジェクト
 */
function createParticle(position={x:0,y:0,z:0}, is_wall=false, velocity={x:0,y:0,z:0}, force={x:0,y:0,z:0}, acceleration={x:0,y:0,z:0}, density=0, pressure=0) {
    return {
        position: {x: position.x, y: position.y, z: position.z},
        velocity: {x: velocity.x, y:velocity.y, z: velocity.z},
        force: {x: force.x, y: force.y, z: force.z},
        acceleration: {x: acceleration.x, y: acceleration.y, z: acceleration.z},
        density: density,
        pressure: pressure,
        is_wall: is_wall
    };
}

/**
 * ベクトルを表すオブジェクトを返す関数
 * @param {number} x - x成分
 * @param {number} y - y成分
 * @param {number} z - z成分
 * @returns {Vector3} - ベクトルを表すオブジェクト
 */
function createVector3(x = 0, y = 0, z = 0) {
    return {x: x, y: y, z: z};
}

/**
 * ベクトルについて a + b を計算する関数
 * @param {Vector3} a - 足すベクトル
 * @param {Vector3} b - 足されるベクトル
 * @returns {Vector3} - ベクトルの和
 */
function addVector3(a, b) {
    return createVector3(a.x + b.x, a.y + b.y, a.z + b.z);
}

/**
 * ベクトルについて a - b を計算する関数
 * @param {Vector3} a - 引かれるベクトル
 * @param {Vector3} b - 引くベクトル
 * @returns {Vector3} - ベクトルの差
 */
function subVector3(a, b) {
    return createVector3(a.x - b.x, a.y - b.y, a.z - b.z);
}

/**
 * ベクトルについて aのn倍 を計算する関数
 * @param {Vector3} a - ベクトル
 * @param {number} n - 実数（スカラー）
 * @returns {Vector3} - 実数倍したベクトル
 */
function multiplyScalarVector3(a, n) {
    return createVector3(n * a.x, n * a.y, n * a.z);
}

/**
 * ベクトルについて a ・ b （内積） を計算する関数
 * @param {Vector3} a - ベクトル
 * @param {Vector3} b - ベクトル
 * @returns {number} - ベクトルの内積
 */
function dotVector3(a, b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}



/**
 * クライアントへ進捗を報告する関数
 * @param {string} text - 報告内容
 */
function reportProgress(text) {
    self.postMessage({type:"progress", content:text});
}